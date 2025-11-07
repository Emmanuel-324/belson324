#!/usr/bin/env python3
"""
particle_stats.py

Extract precipitate geometry from a 2D phase-field Exodus file
for:
  - eta_pv1  (gamma'')
  - eta_pv2  (gamma')
  - eta_pv3  (gamma'')

For each detected precipitate, we record:
- phase_label (gamma_pp or gamma_p)
- precip_id   (ID within that phase map)
- sim_time    (physical simulation time from Exodus: time_whole[-1])
- area_nm2
- centroid_x_nm, centroid_y_nm
- semi_major_nm, semi_minor_nm  (equivalent ellipse from 2nd central moments)
- aspect_ratio = semi_major_nm / semi_minor_nm
- radius_nm    (equivalent-circle radius sqrt(A/pi); useful for gamma')

We follow MacSleyne et al., Acta Materialia 56 (2008) 427–437,
moment-invariants / second-moment-based ellipse extraction.

Output Excel workbook (.xlsx) with four sheets:
  - "gamma_pp_1" : all precipitates from eta_pv1
  - "gamma_p"    : all precipitates from eta_pv2
  - "gamma_pp_3" : all precipitates from eta_pv3
  - "ALL"        : concatenation of all above
"""

import argparse
import numpy as np
import pandas as pd
from scipy.ndimage import label
import meshio
from netCDF4 import Dataset


def compute_region_geometry(x_coords, y_coords, mask_region):
    """
    Given:
      x_coords, y_coords: 2D arrays [nm]
      mask_region: boolean 2D array = True inside precipitate

    Returns dict with:
      area_nm2
      centroid_x_nm, centroid_y_nm
      semi_major_nm, semi_minor_nm, aspect_ratio
      radius_nm
    """

    # Estimate pixel area assuming structured regular grid
    dxs = np.diff(x_coords[0, :])
    dys = np.diff(y_coords[:, 0])
    dx = np.mean(np.abs(dxs))
    dy = np.mean(np.abs(dys))
    dA = dx * dy  # nm^2 per pixel

    idx = np.argwhere(mask_region)
    if idx.size == 0:
        return None

    # Area
    A = idx.shape[0] * dA  # nm^2

    # Pixel coordinates in this precipitate
    xs = x_coords[mask_region]
    ys = y_coords[mask_region]

    # Centroid
    x_c = xs.mean()
    y_c = ys.mean()

    # Central second moments
    x_rel = xs - x_c
    y_rel = ys - y_c

    l20 = np.sum((x_rel**2) * dA)
    l02 = np.sum((y_rel**2) * dA)
    l11 = np.sum((x_rel * y_rel) * dA)

    # Build the moment tensor
    M = np.array([[l20, l11],
                  [l11, l02]], dtype=float)

    # Eigenvalues (principal 2nd moments)
    evals, _ = np.linalg.eigh(M)  # ascending order
    lam_minor, lam_major = evals[0], evals[1]

    # Convert to equivalent ellipse semi-axes a,b using MacSleyne et al.
    # a ~ semi-major axis length, b ~ semi-minor
    if lam_minor <= 0 or lam_major <= 0:
        a = 0.0
        b = 0.0
    else:
        prefac = (4.0 / np.pi) ** 0.25
        a = prefac * ((lam_major**3 / lam_minor) ** 0.125)
        b = prefac * ((lam_minor**3 / lam_major) ** 0.125)
        # enforce ordering
        if b > a:
            a, b = b, a

    aspect_ratio = a / b if (a > 0 and b > 0) else np.nan

    # Equivalent-circle radius (good for gamma' which is closer to equiaxed)
    R_eff = np.sqrt(A / np.pi)

    return dict(
        area_nm2=A,
        centroid_x_nm=x_c,
        centroid_y_nm=y_c,
        semi_major_nm=a,
        semi_minor_nm=b,
        aspect_ratio=aspect_ratio,
        radius_nm=R_eff,
    )


def analyze_phase(x_grid, y_grid, eta_grid, phase_label, sim_time, thr=0.5):
    """
    Segment precipitates for one eta field.
    Returns a list of dict rows for that field only.
    Each row has precip_id starting at 1 for this field.
    """
    if eta_grid is None:
        return []

    # Threshold |eta| > thr to mark precipitate interior
    mask = np.abs(eta_grid) > thr

    # Connected-component labeling
    labeled, nlabels = label(mask)

    rows = []
    for lbl in range(1, nlabels + 1):
        precip_mask = (labeled == lbl)
        geom = compute_region_geometry(x_grid, y_grid, precip_mask)
        if geom is None:
            continue

        row = {
            "phase_label": phase_label,
            "precip_id": int(lbl),
            "sim_time": float(sim_time),
        }
        row.update(geom)
        rows.append(row)

    return rows


def load_exodus_2d_as_grid(exo_path):
    """
    Read:
      - node coordinates
      - eta_pv1, eta_pv2, eta_pv3 at final stored state
      - last simulation time from 'time_whole'

    Assumes:
      - 2D structured mesh (rectangular grid);
      - Exodus file contains multiple times, meshio gives us last state.
    """
    msh = meshio.read(exo_path)

    pts = msh.points  # shape (Npts, 3)
    x_all = pts[:, 0]
    y_all = pts[:, 1]

    xs_unique = np.unique(x_all)
    ys_unique = np.unique(y_all)
    nx = xs_unique.size
    ny = ys_unique.size

    x_grid = np.zeros((ny, nx))
    y_grid = np.zeros((ny, nx))

    # Assemble (ny,nx) coordinate grids by matching coords
    for i, yv in enumerate(ys_unique):
        for j, xv in enumerate(xs_unique):
            hits = np.where(np.isclose(x_all, xv) & np.isclose(y_all, yv))[0]
            if hits.size == 0:
                raise RuntimeError("Could not reconstruct structured grid from Exodus points.")
            idx0 = hits[0]
            x_grid[i, j] = x_all[idx0]
            y_grid[i, j] = y_all[idx0]

    def reshape_field(field_name):
        if field_name not in msh.point_data:
            return None
        arr = msh.point_data[field_name]  # len = Npts
        out = np.zeros((ny, nx))
        for i, yv in enumerate(ys_unique):
            for j, xv in enumerate(xs_unique):
                hits = np.where(np.isclose(x_all, xv) & np.isclose(y_all, yv))[0]
                idx0 = hits[0]
                out[i, j] = arr[idx0]
        return out

    eta1 = reshape_field("eta_pv1")  # gamma'' variant 1
    eta2 = reshape_field("eta_pv2")  # gamma'
    eta3 = reshape_field("eta_pv3")  # gamma'' variant 2

    # Read last simulation time from NetCDF
    try:
        with Dataset(exo_path, "r") as nc:
            times = nc.variables["time_whole"][:]
            sim_time = float(times[-1])
    except Exception as e:
        print(f"Warning: could not read simulation time: {e}")
        sim_time = np.nan

    return x_grid, y_grid, eta1, eta2, eta3, sim_time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--exo", required=True, help="Path to Exodus .e file")
    parser.add_argument("--out", required=True,
                        help="output Excel filename, e.g. particle_stats.xlsx")
    parser.add_argument("--thr", type=float, default=0.5,
                        help="threshold on abs(eta)")
    args = parser.parse_args()

    # Load mesh, fields, physical time
    xg, yg, eta1, eta2, eta3, sim_time = load_exodus_2d_as_grid(args.exo)

    # Analyze each phase separately
    rows_eta1 = analyze_phase(xg, yg, eta1, "gamma_pp_1", sim_time, thr=args.thr)
    rows_eta2 = analyze_phase(xg, yg, eta2, "gamma_p",    sim_time, thr=args.thr)
    rows_eta3 = analyze_phase(xg, yg, eta3, "gamma_pp_3", sim_time, thr=args.thr)

    # Convert each list to its own DataFrame
    cols = [
        "phase_label",
        "precip_id",
        "sim_time",
        "area_nm2",
        "centroid_x_nm",
        "centroid_y_nm",
        "semi_major_nm",
        "semi_minor_nm",
        "aspect_ratio",
        "radius_nm",
    ]

    df_eta1 = pd.DataFrame(rows_eta1, columns=cols)
    df_eta2 = pd.DataFrame(rows_eta2, columns=cols)
    df_eta3 = pd.DataFrame(rows_eta3, columns=cols)

    # ALL sheet = vertically stacked
    df_all = pd.concat([df_eta1, df_eta2, df_eta3], ignore_index=True)

    # Write Excel with 4 sheets:
    #   gamma_pp_1, gamma_p, gamma_pp_3, ALL
    with pd.ExcelWriter(args.out, engine="xlsxwriter") as writer:
        df_eta1.to_excel(writer, sheet_name="gamma_pp_1", index=False)
        df_eta2.to_excel(writer, sheet_name="gamma_p",    index=False)
        df_eta3.to_excel(writer, sheet_name="gamma_pp_3", index=False)
        df_all.to_excel(writer,  sheet_name="ALL",        index=False)

    print(
        f"Wrote {args.out} with "
        f"{len(df_eta1)} (eta_pv1) + {len(df_eta2)} (eta_pv2) + {len(df_eta3)} (eta_pv3) "
        f"= {len(df_all)} total precipitates."
    )
    print(f"Simulation time recorded: {sim_time}")


if __name__ == "__main__":
    main()
