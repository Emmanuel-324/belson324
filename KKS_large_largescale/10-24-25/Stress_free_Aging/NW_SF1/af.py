#!/usr/bin/env python3
import pandas as pd
import numpy as np
from pathlib import Path

# --- Settings ---
domain_length_nm = 350.0
domain_area_nm2 = domain_length_nm ** 2

# --- Input file path ---
file_path = Path("/home/emmanuel324/projects/belson324/KKS_large_largescale/10-24-25/Stress_free_Aging/NW_SF1/AR_NW_SF1.xlsx")

# --- Output file path (same folder, different name) ---
output_path = file_path.with_name(file_path.stem + "_areafraction.xlsx")

# --- Read Excel sheet (γʺ data) ---
df = pd.read_excel(file_path, sheet_name="gamma_pp_ALL")

# --- Compute area of each precipitate (ellipse area = πab) ---
df['area_nm2'] = np.pi * df['mean_major_nm'] * df['mean_minor_nm']

# --- Group by step or time (whichever exists) ---
if 'step_index' in df.columns:
    grouped = df.groupby('step_index')
elif 'time_s' in df.columns:
    grouped = df.groupby('time_s')
else:
    grouped = [(0, df)]  # single snapshot

# --- Compute area fraction for each group ---
area_fraction = []
for step, g in grouped:
    total_precip_area = g['area_nm2'].sum()
    frac_percent = 100 * total_precip_area / domain_area_nm2
    area_fraction.append({
        'step_index': step,
        'area_fraction_percent': frac_percent
    })

# --- Convert to DataFrame and save ---
area_df = pd.DataFrame(area_fraction)
area_df.to_excel(output_path, index=False)

print(f"✅ Area fraction results saved to:\n{output_path}")
print(area_df)
