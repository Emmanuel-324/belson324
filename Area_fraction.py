import pandas as pd
import numpy as np
import re
from pathlib import Path
from typing import List, Dict, Any

# ========= Settings =========
file_path = Path(r"/home/emmanuel324/projects/belson324/KKS_large_largescale/9-12-25/ES_04/ES_04_table.csv")  # your CSV
# Output paths: if None, will write next to the input file with suffixes
output_csv: Path | None = None
# If True, include the time_hours column in the output alongside area fractions
include_time_in_output: bool = True
output_xlsx: Path | None = Path(r"/home/emmanuel324/projects/belson324/KKS_large_largescale/9-12-25/ES_04/Area_fraction_ES04.xlsx")
domain_length_x_nm = 350
domain_length_y_nm = 350

# Auto-detect all eta_pv* columns by default; or set explicitly e.g. ['eta_pv1', 'eta_pv2', 'eta_pv3']
precip_columns: List[str] = None

# Map precip columns to output names for percentage area fraction (D1/D2/D3)
NAME_MAP = {
    'eta_pv1': 'D1',
    'eta_pv2': 'D2',
    'eta_pv3': 'D3',
}

# Per-variant downsampling config (exact row counts for Excel)
DOWNSAMPLE_CONFIG: Dict[str, Dict[str, Any]] = {
    'eta_pv1': {'mode': 'n_points', 'n_points': 20},
    'eta_pv2': {'mode': 'n_points', 'n_points': 20},
    'eta_pv3': {'mode': 'n_points', 'n_points': 20},
}

# ========= Helpers =========
def even_spaced_indices_exact(n_total: int, n_points: int) -> np.ndarray:
    """
    Return exactly min(n_points, n_total) strictly increasing indices from [0, n_total-1],
    spread as evenly as possible and including the endpoints when feasible.
    """
    k = int(max(1, n_points))
    if n_total <= k:
        return np.arange(n_total, dtype=int)

    # initial evenly-spaced float positions (0 .. n_total-1)
    pos = np.linspace(0, n_total - 1, k)
    idx = np.round(pos).astype(int)

    # enforce strict monotonic increase
    idx[0] = max(0, idx[0])
    for i in range(1, k):
        if idx[i] <= idx[i - 1]:
            idx[i] = idx[i - 1] + 1

    # clamp to bounds [0, n_total-1]
    if idx[-1] > n_total - 1:
        # shift backward to fit into the range without duplicates
        # make it an arithmetic sequence ending at n_total-1
        idx = np.arange(n_total - k, n_total, dtype=int)

    return idx

def downsample_n_points(df: pd.DataFrame, n_points: int) -> pd.DataFrame:
    """
    Downsample dataframe rows to exactly min(n_points, len(df)) rows using even spacing.
    Preserves all columns.
    """
    n = len(df)
    indices = even_spaced_indices_exact(n, n_points)
    return df.iloc[indices].copy()

# ========= Load =========
def robust_read_csv(path: Path) -> pd.DataFrame:
    """Try to read CSV with the fast engine first, fall back to python engine with
    relaxed bad-line handling if a ParserError occurs. Print diagnostics to help
    locate malformed lines.
    """
    try:
        return pd.read_csv(path)
    except pd.errors.ParserError as e:
        print(f"ParserError reading {path}: {e}")
        print("Falling back to python engine with on_bad_lines='warn' (will skip/break bad rows).")
        try:
            df = pd.read_csv(path, engine='python', on_bad_lines='warn')
            return df
        except Exception as e2:
            print(f"Fallback read also failed: {e2}")
            # Try to show the offending line numbers by scanning the file for inconsistent field counts
            try:
                with path.open('r', encoding='utf-8', errors='replace') as f:
                    first = None
                    problematic = []
                    for i, line in enumerate(f, start=1):
                        if not line.strip():
                            continue
                        fields = line.rstrip('\n').split(',')
                        if first is None:
                            first = len(fields)
                        elif len(fields) != first:
                            problematic.append((i, len(fields), line.strip()))
                            if len(problematic) >= 10:
                                break
                if problematic:
                    print('Detected lines with unexpected field counts (first 10):')
                    for ln, cnt, text in problematic:
                        print(f"  Line {ln}: {cnt} fields -> {text[:200]}")
                else:
                    print('No inconsistent line lengths found when scanning; file may be using a different delimiter or quoting scheme.')
            except Exception as e3:
                print('While scanning file for diagnostics an error occurred:', e3)
            raise


df = robust_read_csv(file_path)
# Clean column names: strip surrounding whitespace and collapse internal whitespace to single underscore
df.columns = df.columns.str.strip().str.replace(r"\s+", "_", regex=True)

# ========= Pick precipitate columns (robust) =========
# Build a normalized map of column name -> normalized identifier (lowercase, stripped, alnum + underscore)
col_norm_map = {c: re.sub(r'[^a-z0-9_]', '', c.strip().lower()) for c in df.columns}

if precip_columns is None:
    # Detect columns whose normalized name starts with 'eta' (covers 'eta_pv1', 'eta1', etc.)
    precip_columns = [orig for orig, norm in col_norm_map.items() if norm.startswith('eta')]

if not precip_columns:
    # helpful diagnostic
    available = ", ".join(list(df.columns))
    raise ValueError(
        "No precipitate columns found (looking for names starting with 'eta' or 'eta_pv').\n"
        f"Available columns: {available}\n"
        "If your eta columns have leading/trailing spaces or unusual characters, please clean them or set `precip_columns` explicitly."
    )

# We'll use normalized keys for NAME_MAP lookups. Build helper to get normalized key for an original column name.
def _norm(col_name: str) -> str:
    return col_norm_map.get(col_name, re.sub(r'[^a-z0-9_]', '', col_name.strip().lower()))

# ========= Time: seconds -> hours (add, don't replace) =========
if 'time' in df.columns:
    df['time_hours'] = df['time'] / 3600.0
else:
    df['time_hours'] = df.index.astype(float)

# ========= Compute per-variant area fractions (percent) into D1/D2/D3 =========
total_material_area = float(domain_length_x_nm) * float(domain_length_y_nm)

# Remove any old D-columns we are about to write, then recompute fresh
for col in precip_columns:
    norm = _norm(col)
    # prefer normalized NAME_MAP key, fall back to a safe generated name without spaces
    mapped = NAME_MAP.get(norm) or NAME_MAP.get(col)
    if mapped:
        dname = mapped
    else:
        safe_col = re.sub(r'[^0-9a-zA-Z_]', '_', col.strip())
        dname = f"D_{safe_col}"

    if dname in df.columns:
        del df[dname]

    # use the original column name to read values, but write the standardized D-name
    df[dname] = (df[col] / total_material_area) * 100.0

# ========= Save full dataframe to a NEW CSV/XLSX (do not overwrite input by default) =========
if output_csv is None:
    output_csv = file_path.with_name(file_path.stem + "_with_areafractions.csv")

# Build a minimal output dataframe containing only the area fraction columns (and optionally time_hours)
# Derive area columns from the dataframe itself: any column whose name starts with 'D' and was just created above
area_cols = [c for c in df.columns if str(c).startswith('D')]
output_columns = []
if include_time_in_output and 'time_hours' in df.columns:
    output_columns.append('time_hours')
output_columns.extend(area_cols)

out_df = df.loc[:, output_columns].copy()

out_df.to_csv(output_csv, index=False)
print("Area-fraction CSV written to:", output_csv)

# Optionally also write an Excel version containing only the area fractions
if output_xlsx is None:
    output_xlsx = file_path.with_name(file_path.stem + "_with_areafractions.xlsx")

try:
    out_df.to_excel(output_xlsx, index=False)
    print("Area-fraction Excel written to:", output_xlsx)
except Exception as e:
    # If openpyxl / engine not available, warn but continue (downsampled Excel is still saved below)
    print("Warning: failed to write area-fraction Excel file:", e)

# ========= Build per-variant downsampled frames (EXACT row counts) =========
downsampled_frames: Dict[str, pd.DataFrame] = {}

for col in precip_columns:
    cfg = DOWNSAMPLE_CONFIG.get(col, {'mode': 'n_points', 'n_points': 20})
    mode = (cfg.get('mode') or 'n_points').lower()
    dname = NAME_MAP.get(col, f"D_{col}")  # 'D1', 'D2', 'D3', ...

    if mode != 'n_points':
        raise ValueError(f"For this script, only 'n_points' mode is enforced. Got mode='{mode}' for {col}.")

    n_points = int(max(1, cfg.get('n_points', 20)))
    ds_df = downsample_n_points(df, n_points)
    # keep ALL columns (“the rest untouched”), but only the selected rows
    downsampled_frames[dname] = ds_df

    # sanity check: exact row count
    assert len(ds_df) == min(n_points, len(df)), \
        f"{dname}: expected {min(n_points, len(df))} rows, got {len(ds_df)}"

    print(f"{dname}: downsampled rows saved = {len(ds_df)}")

# ========= Save ONLY the DOWNSAMPLED results to a single Excel file (sheets D1/D2/D3) =========
out_xlsx = file_path.with_name(file_path.stem + "_downsampled.xlsx")
with pd.ExcelWriter(out_xlsx) as writer:
    for dname, ds_df in downsampled_frames.items():
        ds_df.to_excel(writer, sheet_name=dname, index=False)

print("Downsampled Excel saved to:", out_xlsx)
print("Sheets:", list(downsampled_frames.keys()))
