from pathlib import Path
import argparse
import pandas as pd

# ----- Config you can tweak if needed -----
R_MIN_NM = 2.0   # γ' (eta_pv2) radius floor for PASSED tables

# ----- Args -----
ap = argparse.ArgumentParser()
ap.add_argument("--csv", required=True, type=Path, help="Full path to aspect_ratios_by_precip.csv")
ap.add_argument("--out", type=Path, default=None, help="Optional output .xlsx (default: alongside CSV)")
ap.add_argument("--start_step", type=int, default=0, help="Ignore rows with step_index < START (e.g., 85)")
args = ap.parse_args()

csv_path = args.csv
out_xlsx = args.out if args.out else csv_path.with_name("summary_all_vs_passed.xlsx")

# ----- Load -----
df = pd.read_csv(csv_path)
print(f"Loaded: {csv_path}")

# Required columns
for c in ["phase","variant","step_index","time_value"]:
    if c not in df.columns:
        raise ValueError(f"Missing column '{c}' in input CSV.")
if "passes_Lmin_nm" not in df.columns:
    df["passes_Lmin_nm"] = True  # fallback if older CSVs lack the flag

# ----- Apply start-step cutoff only -----
df = df[df["step_index"] >= args.start_step].copy()
print(f"Filtering: keeping step_index >= {args.start_step}. Remaining rows: {len(df)}")

# ----- Helpers (no per-step count gating) -----
def summarize_gamma_pp(table: pd.DataFrame) -> pd.DataFrame:
    if table.empty:
        return pd.DataFrame(columns=[
            "step_index","time_s","n",
            "mean_minor_nm","median_minor_nm",
            "mean_major_nm","median_major_nm",
            "mean_AR","median_AR"
        ])
    g = (table.groupby("step_index")
               .agg(time_s=("time_value","median"),
                    n=("minor_nm","size"),
                    mean_minor_nm=("minor_nm","mean"),
                    median_minor_nm=("minor_nm","median"),
                    mean_major_nm=("major_nm","mean"),
                    median_major_nm=("major_nm","median"),
                    mean_AR=("aspect_ratio","mean"),
                    median_AR=("aspect_ratio","median"))
               .reset_index())
    return g.sort_values("step_index", kind="stable").reset_index(drop=True)

def summarize_gamma_p(table: pd.DataFrame) -> pd.DataFrame:
    if table.empty:
        return pd.DataFrame(columns=["step_index","time_s","n","mean_radius_nm","median_radius_nm"])
    g = (table.groupby("step_index")
               .agg(time_s=("time_value","median"),
                    n=("radius_nm","size"),
                    mean_radius_nm=("radius_nm","mean"),
                    median_radius_nm=("radius_nm","median"))
               .reset_index())
    return g.sort_values("step_index", kind="stable").reset_index(drop=True)

# ----- Build four summaries (ALL vs PASSED) -----
# γʺ = gamma_pp (eta_pv1 + eta_pv3)
gpp_all  = df[df["phase"] == "gamma_pp"].copy()
gpp_pass = df[(df["phase"] == "gamma_pp") & (df["passes_Lmin_nm"] == True)].copy()
gamma_pp_ALL    = summarize_gamma_pp(gpp_all)
gamma_pp_PASSED = summarize_gamma_pp(gpp_pass)

# γ' = gamma_p (eta_pv2)
gp_all  = df[df["phase"] == "gamma_p"].copy()
gp_pass = df[(df["phase"] == "gamma_p") & (df["radius_nm"] >= R_MIN_NM)].copy()
gamma_p_ALL     = summarize_gamma_p(gp_all)
gamma_p_PASSED  = summarize_gamma_p(gp_pass)

# ----- Save -----
with pd.ExcelWriter(out_xlsx) as xw:
    gamma_pp_ALL.to_excel(xw, sheet_name="gamma_pp_ALL", index=False)
    gamma_pp_PASSED.to_excel(xw, sheet_name="gamma_pp_PASSED", index=False)
    gamma_p_ALL.to_excel(xw, sheet_name="gamma_p_ALL", index=False)
    gamma_p_PASSED.to_excel(xw, sheet_name="gamma_p_PASSED", index=False)

print("Saved:", out_xlsx)
print({
    "gamma_pp_ALL_first_step": None if gamma_pp_ALL.empty else int(gamma_pp_ALL["step_index"].min()),
    "gamma_pp_PASSED_first_step": None if gamma_pp_PASSED.empty else int(gamma_pp_PASSED["step_index"].min()),
    "gamma_p_ALL_first_step": None if gamma_p_ALL.empty else int(gamma_p_ALL["step_index"].min()),
    "gamma_p_PASSED_first_step": None if gamma_p_PASSED.empty else int(gamma_p_PASSED["step_index"].min()),
})
