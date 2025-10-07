#!/usr/bin/env python3
# Compute 2D area fractions (%) per time step from an Exodus .e file for precipitate order parameters.
# This script assumes the phase indicators are order parameters (eta_pv1, eta_pv2, eta_pv3) and uses abs() for fractions since etas can be negative.
# Requires: pip install netCDF4

from pathlib import Path
from typing import Dict, List, Tuple
import argparse
import numpy as np
from netCDF4 import Dataset
import csv

def decode_names(char_arr) -> List[str]:
    if char_arr.ndim != 2:
        return []
    out = []
    # Handle MaskedArray by getting the underlying data
    char_data = np.ma.getdata(char_arr) if np.ma.is_masked(char_arr) else char_arr
    for i in range(char_data.shape[0]):
        row = char_data[i]
        if row.dtype.kind == 'S':  # Bytes/string dtype (e.g., |S1 or |S8)
            # Join bytes elements into a single bytes object, then decode
            byte_row = b''.join(row.astype('|S1'))
            s = byte_row.decode('utf-8').strip('\x00')  # Strip null bytes too
        else:  # Assume integer dtype (ASCII codes)
            # Convert each element to int and then to chr, skipping any invalid
            valid_chars = [c for c in row if 0 <= c < 256]  # Filter valid byte range
            s = ''.join(chr(int(c)) for c in valid_chars).strip('\x00')
        # Only append non-empty strings
        if s:
            out.append(s)
    return out

def polygon_area_xy(x: np.ndarray, y: np.ndarray) -> float:
    # Shoelace formula
    return 0.5 * np.abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))

def collect_mesh(ds: Dataset):
    if "coordx" not in ds.variables or "coordy" not in ds.variables:
        raise RuntimeError("Expected 'coordx' and 'coordy' for a 2D Exodus mesh.")
    X = ds.variables["coordx"][:]
    Y = ds.variables["coordy"][:]

    if "num_el_blk" not in ds.dimensions:
        raise RuntimeError("No element blocks found in Exodus file.")
    nblk = ds.dimensions["num_el_blk"].size

    elems_per_block: List[np.ndarray] = []
    for b in range(nblk):
        blk_id = b + 1
        name = f"connect{blk_id}"
        if name not in ds.variables:
            raise RuntimeError(f"Missing connectivity var '{name}'.")
        conn = ds.variables[name][:].astype(int)  # 1-based node IDs
        elems_per_block.append(conn)
    return X, Y, elems_per_block

def areas_per_element(X, Y, elems_per_block) -> np.ndarray:
    areas = []
    for conn in elems_per_block:
        for nodes in conn:
            ids0 = nodes - 1
            x = X[ids0]; y = Y[ids0]
            areas.append(polygon_area_xy(x, y))
    return np.asarray(areas)

def read_times(ds: Dataset) -> np.ndarray:
    if "time_whole" in ds.variables:
        return ds.variables["time_whole"][:]
    if "time" in ds.variables:
        return ds.variables["time"][:]
    raise RuntimeError("No 'time_whole' or 'time' variable found in Exodus file.")

def name_index_map(char_var_name: str, ds: Dataset) -> Dict[str, int]:
    if char_var_name not in ds.variables:
        return {}
    names = decode_names(ds.variables[char_var_name][:])
    return {n.lower(): i for i, n in enumerate(names)}

def elem_var_handle(ds: Dataset, var_idx0: int, blk_id: int):
    name = f"vals_elem_var{var_idx0+1}eb{blk_id}"
    return ds.variables.get(name, None)

def get_element_values_for_var(ds: Dataset,
                               varname: str,
                               elems_per_block,
                               name_to_idx_elem,
                               name_to_idx_nod,
                               t_idx0: int) -> np.ndarray:
    key = varname.lower()

    # Prefer element vars
    if key in name_to_idx_elem:
        iv0 = name_to_idx_elem[key]
        blocks = []
        for b, conn in enumerate(elems_per_block):
            vh = elem_var_handle(ds, iv0, b + 1)
            if vh is None:
                raise RuntimeError(f"Element data not found for '{varname}' in block {b+1}")
            blocks.append(vh[t_idx0, :])
        return np.concatenate(blocks, axis=0)

    # Fallback: nodal var -> average to element centers
    if key in name_to_idx_nod:
        iv0 = name_to_idx_nod[key]
        vnod = ds.variables[f"vals_nod_var{iv0+1}"][t_idx0, :]  # [num_nodes]
        vals = []
        for conn in elems_per_block:
            for nodes in conn:
                ids0 = nodes - 1
                vals.append(float(np.mean(vnod[ids0])))
        return np.asarray(vals)

    raise RuntimeError(f"Variable '{varname}' not found as element or nodal variable.")

def compute_area_fractions_percent(exo_file: Path, phase_vars: List[str]):
    with Dataset(exo_file, "r") as ds:
        times = read_times(ds)
        X, Y, elems_per_block = collect_mesh(ds)
        areas = areas_per_element(X, Y, elems_per_block)
        A_total = float(np.sum(areas))

        name_to_idx_elem = name_index_map("name_elem_var", ds)
        name_to_idx_nod  = name_index_map("name_nod_var", ds)

        out: Dict[str, np.ndarray] = {}
        for pv in phase_vars:
            arr = np.zeros_like(times, dtype=float)
            for it in range(len(times)):
                h_e = get_element_values_for_var(ds, pv, elems_per_block,
                                                 name_to_idx_elem, name_to_idx_nod, it)
                arr[it] = 100.0 * float(np.dot(np.abs(h_e), areas) / A_total)  # Use abs() for order parameters, percent
            out[pv] = arr
        return times, out

def save_csv(out_path: Path, times: np.ndarray, fracs: Dict[str, np.ndarray], header_names: List[str]) -> Path:
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time"] + header_names)
        n = len(times)
        for i in range(n):
            row = [times[i]] + [fracs[k][i] for k in header_names]
            w.writerow(row)
    return out_path

def main():
    ap = argparse.ArgumentParser(description="Compute 2D area fractions (%) from Exodus .e file for precipitate order parameters (using abs for negatives).")
    ap.add_argument("--exo", required=True, type=Path, help="Path to Exodus .e file (e.g., ES_01a_out.e)")
    ap.add_argument("--vars", nargs="+", default=["eta_pv1", "eta_pv2", "eta_pv3"], help="Order parameter variable names (default: eta_pv1 eta_pv2 eta_pv3)")
    ap.add_argument("--out", type=Path, default=None, help="Optional output CSV path (default: next to .e)")
    args = ap.parse_args()

    exo = args.exo
    if not exo.exists():
        raise FileNotFoundError(f"Exodus file not found: {exo}")

    phase_vars = args.vars
    times, fracs = compute_area_fractions_percent(exo, phase_vars)

    out_csv = args.out if args.out is not None else exo.with_name(exo.stem + "_area_fractions_percent.csv")
    # Preserve given order of --vars in the CSV columns
    save_csv(out_csv, times, fracs, phase_vars)

    # Quick sanity: sum ~ 100% (may not include matrix eta_m)
    sums = np.zeros_like(times, dtype=float)
    for k in phase_vars:
        sums += fracs[k]
    print(f"Wrote: {out_csv}")
    print("Columns:", ["time"] + phase_vars)
    print("Mean(sum of % area fractions) ≈", float(np.mean(sums)))

if __name__ == "__main__":
    main()