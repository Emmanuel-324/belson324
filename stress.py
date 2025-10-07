import pandas as pd
from pathlib import Path
import numpy as np

# ---- SET THIS ----
file_path = Path(r"C:/Users/adugy/Downloads/9-19-25/ES_03/ES_03_table.csv")  # or .xlsx

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

def write_any(df: pd.DataFrame, path: Path):
    if path.suffix.lower() in [".xlsx", ".xls"]:
        # overwrite same Excel file (single sheet)
        with pd.ExcelWriter(path, engine="openpyxl", mode="w") as xw:
            df.to_excel(xw, index=False)
    else:
        df.to_csv(path, index=False)

df = read_any(file_path)

required = ["num_pv1","num_pv2","num_pv3","num_m",
            "den_pv1","den_pv2","den_pv3","den_m"]
missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns: {missing}")

# safe division helper: returns NaN when denominator = 0
def safe_div(n, d):
    n = pd.to_numeric(n, errors="coerce")
    d = pd.to_numeric(d, errors="coerce")
    return np.where(d != 0, n / d, np.nan)

df["avg_pv1_stress_xx"] = safe_div(df["num_pv1"], df["den_pv1"])
df["avg_pv2_stress_xx"] = safe_div(df["num_pv2"], df["den_pv2"])
df["avg_pv3_stress_xx"] = safe_div(df["num_pv3"], df["den_pv3"])
df["avg_m_stress_xx"]   = safe_div(df["num_m"],   df["den_m"])

write_any(df, file_path)
print("Per-phase average stress (σ_xx) added:",
      "[avg_pv1_stress_xx, avg_pv2_stress_xx, avg_pv3_stress_xx, avg_m_stress_xx]")
print(f"Saved to: {file_path}")

