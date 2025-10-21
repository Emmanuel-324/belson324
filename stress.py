import pandas as pd
from pathlib import Path
import numpy as np

# ---- SET THIS ----
file_path = Path(r"/home/emmanuel324/projects/belson324/KKS_large_largescale/10-13-25/EigenA1_a/EigenA1_a_table.csv")  # or .xlsx

# ------------------

def read_any(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in [".xlsx", ".xls"]:
        return pd.read_excel(path)
    else:
        # default to CSV (handles comma or tab if needed)
        try:
            return pd.read_csv(path)
        except Exception:
            return pd.read_csv(path, sep="\t")

df = read_any(file_path)

required = ["num_vonmises_pv1","num_vonmises_pv2","num_vonmises_pv3","num_vonmises_m",
            "den_pv1","den_pv2","den_pv3","den_m"]
missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns: {missing}")

# safe division helper: returns NaN when denominator = 0
def safe_div(n, d):
    n = pd.to_numeric(n, errors="coerce")
    d = pd.to_numeric(d, errors="coerce")
    return np.where(d != 0, n / d, np.nan)

# Calculate von Mises stress for each phase
df["pv1_vonmises_stress"] = safe_div(df["num_vonmises_pv1"], df["den_pv1"])
df["pv2_vonmises_stress"] = safe_div(df["num_vonmises_pv2"], df["den_pv2"])
df["pv3_vonmises_stress"] = safe_div(df["num_vonmises_pv3"], df["den_pv3"])
df["m_vonmises_stress"]   = safe_div(df["num_vonmises_m"],   df["den_m"])

# Create the output path for the new Excel file in the same directory
output_dir = file_path.parent
output_path = output_dir / "stress.xlsx"

# Write the results to a new Excel file
with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
    # Write the full dataset with new columns
    df.to_excel(writer, sheet_name="stress_analysis", index=False)
    
    # Create a summary sheet with just the stress columns
    stress_cols = ['time'] + [col for col in df.columns if 'vonmises_stress' in col]
    df[stress_cols].to_excel(writer, sheet_name="stress_summary", index=False)

print("Per-phase von Mises stress (σ_xx) calculated:")
print("  - pv1_vonmises_stress: von Mises stress in γ'' variant 1")
print("  - pv2_vonmises_stress: von Mises stress in γ'' variant 2") 
print("  - pv3_vonmises_stress: von Mises stress in γ'' variant 3")
print("  - m_vonmises_stress: von Mises stress in matrix (γ)")
print(f"\nResults saved to: {output_path}")
print(f"  - Sheet 'stress_analysis': Full dataset with stress columns")
print(f"  - Sheet 'stress_summary': Time vs stress columns only")