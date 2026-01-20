#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Aspect ratio / radius summaries with exclusive phase segmentation
for a single Exodus file (e.g. CS_2_out.e) containing all time steps.

Mapping of phase fields (from Exodus):
  eta_gp   = γ'      (roughly spherical → radius statistics)
  eta_gpp1 = γ'' v1  (plate-like → aspect ratio statistics)
  eta_gpp2 = γ'' v2  (plate-like → aspect ratio statistics)

The script:
  - Reads the Exodus file (all time steps)
  - Resamples the three eta fields onto a uniform 2D grid
  - Builds exclusive masks (no overlap) for γ', γ''
  - Measures plate aspect ratios for γ'' (both variants combined)
  - Measures equivalent radius for γ'
  - Writes per-time summaries to Excel (or CSV fallback).
"""

from pathlib import Path
import re
import argparse
import numpy as np
import pandas as pd
import pyvista as pv
from skimage import measure, morphology

# ===================== CLI =====================
ap = argparse.ArgumentParser(
    description="Aspect ratio / radius summaries with exclusive phase segmentation."
)
ap.add_argument(
    "--root", type=Path, required=True,
    help="Directory containing the Exodus file (e.g. CS_2_out.e)"
)
ap.add_argument(
    "--master_glob", type=str, default="*.e",
    help="Glob for main .e file (default: *.e, e.g. CS_2_out.e)"
)
ap.add_argument(
    "--start_step", type=int, default=0,
    help="Ignore rows with step_index < START (e.g., 0, 50, ...)"
)
ap.add_argument(
    "--file_stride", type=int, default=1,
    help="(Legacy) If companion -s### files exist, process every k-th file"
)
ap.add_argument(
    "--max_files", type=int, default=None,
    help="(Legacy) Limit number of companion files"
)
ap.add_argument(
    "--out", type=Path, default=None,
    help="Output .xlsx (or CSV fallback) path; default derived from Exodus file name"
)

# physics / segmentation
ap.add_argument(
    "--L_MIN_NM", type=float, default=3.0,
    help="γʺ: minor-axis pass floor (nm)"
)
ap.add_argument(
    "--R_MIN_NM", type=float, default=2.0,
    help="γʹ: radius pass floor (nm)"
)
ap.add_argument(
    "--tau_gp", type=float, default=0.0,
    help="threshold for eta_gp (γ')"
)
ap.add_argument(
    "--tau_gpp1", type=float, default=0.0,
    help="threshold for eta_gpp1 (γʺ v1)"
)
ap.add_argument(
    "--tau_gpp2", type=float, default=0.0,
    help="threshold for eta_gpp2 (γʺ v2)"
)
ap.add_argument(
    "--open_frac", type=float, default=0.25,
    help="opening radius as a fraction of L_MIN_NM (0.25 = mild)"
)

# units / sampling
ap.add_argument(
    "--unit_to_nm", type=float, default=1.0,
    help="Multiply coordinates by this to get nm (1.0 if already nm)"
)
ap.add_argument(
    "--auto_res", action="store_true",
    help="Auto-pick NX/NY from bounds & L_MIN_NM (~6 px across L_MIN)"
)
ap.add_argument(
    "--nx", type=int, default=200,
    help="Grid NX if not using --auto_res (default: 200)"
)
ap.add_argument(
    "--ny", type=int, default=200,
    help="Grid NY if not using --auto_res (default: 200)"
)

args = ap.parse_args()

root = args.root

# ===================== PyVista compat =====================
UniformGrid = getattr(pv, "UniformGrid", pv.ImageData)


# ===================== Helpers: Exodus =====================
def parse_step_from_name(name: str) -> int:
    """Parse step index from *.e-s### filename (legacy)."""
    m = re.search(r"\.(?:e|exo)-s(\d+)$", name, flags=re.I)
    return int(m.group(1)) if m else 0


def read_time_value_any(ds) -> float | None:
    """Try to read a time value stored in FieldData on dataset or blocks."""
    def grab(fd):
        for k in ("TimeValue", "time", "TIME", "Time"):
            if k in fd:
                v = fd[k]
                try:
                    return float(v[0]) if hasattr(v, "__len__") else float(v)
                except Exception:
                    pass
        return None

    if hasattr(ds, "field_data"):
        t = grab(ds.field_data)
        if t is not None:
            return t

    if isinstance(ds, pv.MultiBlock):
        for ch in ds:
            if ch is not None and hasattr(ch, "field_data"):
                t = grab(ch.field_data)
                if t is not None:
                    return t

    return None


def find_variant_in_dataset(dset, accepted_names):
    """Find one of accepted_names in point_data or cell_data, case-insensitive."""
    pkeys = list(dset.point_data.keys())
    lp = {k.lower(): k for k in pkeys}
    ckeys = list(dset.cell_data.keys())
    lc = {k.lower(): k for k in ckeys}

    for nm in accepted_names:
        if nm in lp:
            return lp[nm], "point"
    for nm in accepted_names:
        if nm in lc:
            return lc[nm], "cell"
    return None, None


def detect_eta_arrays_any(ds):
    """
    Detect η-fields in the Exodus dataset.

    Internal tags in this script:
      eta_pv1 → γʺ (variant 1)
      eta_pv2 → γʹ
      eta_pv3 → γʺ (variant 2)

    Mapping to YOUR Exodus variable names:
      eta_pv1: look for 'eta_gpp1'
      eta_pv2: look for 'eta_gp'
      eta_pv3: look for 'eta_gpp2'
    """
    want = {
        "eta_pv1": ["eta_gpp1", "eta1", "eta_pv1"],
        "eta_pv2": ["eta_gp", "eta2", "eta_pv2"],
        "eta_pv3": ["eta_gpp2", "eta3", "eta_pv3"],
    }
    out = {k: None for k in want}

    def maybe_set(tag, val, blk):
        if val and not out[tag]:
            out[tag] = {"name": val[0], "location": val[1], "block": blk}

    if isinstance(ds, pv.MultiBlock):
        for bi, ch in enumerate(ds):
            if ch is None:
                continue
            if isinstance(ch, pv.MultiBlock):
                for bj, gch in enumerate(ch):
                    if gch is None:
                        continue
                    for tag, alts in want.items():
                        maybe_set(tag, find_variant_in_dataset(gch, alts), (bi, bj))
            else:
                for tag, alts in want.items():
                    maybe_set(tag, find_variant_in_dataset(ch, alts), (bi,))
    else:
        for tag, alts in want.items():
            nm, loc = find_variant_in_dataset(ds, alts)
            if nm:
                out[tag] = {"name": nm, "location": loc, "block": ()}

    return out


def get_block(ds, idx_tuple):
    """Navigate a MultiBlock with an index tuple; empty tuple returns ds."""
    if not idx_tuple:
        return ds
    cur = ds
    for k in idx_tuple:
        cur = cur[k]
    return cur


# ===================== Helpers: sampling & masks =====================
def choose_resolution_from_bounds(bounds, unit_to_nm: float, L_min_nm: float,
                                  min_dim=128, max_dim=2048):
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    xmin *= unit_to_nm
    xmax *= unit_to_nm
    ymin *= unit_to_nm
    ymax *= unit_to_nm
    Lx = max(xmax - xmin, 1e-9)
    Ly = max(ymax - ymin, 1e-9)

    # ~6 pixels across the smallest plate thickness scale
    p_target = max(L_min_nm / 6.0, 1e-6)
    NX = int(np.clip(round(Lx / p_target) + 1, min_dim, max_dim))
    NY = int(np.clip(round(Ly / p_target) + 1, min_dim, max_dim))
    return NX, NY


def resample_to_uniform(block, array_name, location, nx, ny, unit_to_nm):
    """Resample one array to a uniform 2D grid."""
    if location == "cell":
        block = block.cell_data_to_point_data()

    try:
        surf = block.extract_surface()
        work = surf if surf.n_points > 0 else block
    except Exception:
        work = block

    xmin, xmax, ymin, ymax, zmin, zmax = work.bounds
    xmin *= unit_to_nm
    xmax *= unit_to_nm
    ymin *= unit_to_nm
    ymax *= unit_to_nm
    zmin *= unit_to_nm
    zmax *= unit_to_nm

    z0 = 0.5 * (zmin + zmax)
    dx = (xmax - xmin) / max(nx - 1, 1)
    dy = (ymax - ymin) / max(ny - 1, 1)

    grid = UniformGrid()
    grid.dimensions = (nx, ny, 1)
    grid.origin = (xmin, ymin, z0)
    grid.spacing = (dx, dy, max(zmax - zmin, 1.0))

    smp = grid.sample(work)
    if array_name not in smp.point_data:
        raise KeyError(f"Array '{array_name}' missing after sampling.")
    vals = np.asarray(smp.point_data[array_name])
    img = vals.reshape((ny, nx))
    return img, dx, dy, (xmin, xmax, ymin, ymax)


def px_from_nm(L_nm, px_nm_x, px_nm_y):
    return L_nm / max(np.sqrt(px_nm_x * px_nm_y), 1e-12)


def min_area_px_from_minor(L_nm, px_nm_x, px_nm_y):
    mpx = px_from_nm(L_nm, px_nm_x, px_nm_y)
    return int(np.ceil((np.pi / 4.0) * mpx**2))


def exclusive_phase_masks(eta1_img, eta2_img, eta3_img,
                          px_nm_x, px_nm_y,
                          tau1, tau2, tau3,
                          L_MIN_NM, open_frac):
    """
    Build mutually exclusive masks:
      m1: γʺ from eta_pv1 (your eta_gpp1)
      m2: γ'  from eta_pv2 (your eta_gp)
      m3: γʺ from eta_pv3 (your eta_gpp2)
    """
    e1 = np.asarray(eta1_img, float)
    e2 = np.asarray(eta2_img, float)
    e3 = np.asarray(eta3_img, float)

    # Dominance (argmax)
    dom1 = (e1 >= e2) & (e1 >= e3)
    dom2 = (e2 >= e1) & (e2 >= e3)
    dom3 = (e3 >= e1) & (e3 >= e2)

    m1 = (e1 >= tau1) & dom1
    m2 = (e2 >= tau2) & dom2
    m3 = (e3 >= tau3) & dom3

    # Light cleanup: remove tiny specks for γʺ
    numeric_floor_nm = max(1.0, 0.5 * L_MIN_NM)
    min_area_px = min_area_px_from_minor(numeric_floor_nm, px_nm_x, px_nm_y)
    if min_area_px > 0:
        m1 = morphology.remove_small_objects(m1, min_size=min_area_px)
        m3 = morphology.remove_small_objects(m3, min_size=min_area_px)
        # For γ' we can keep very small features; they are usually less problematic.

    # Morphological opening to break skinny bridges
    open_px = max(1, int(px_from_nm(L_MIN_NM, px_nm_x, px_nm_y) * open_frac))
    if open_px > 0:
        se = morphology.disk(open_px)
        m1 = morphology.binary_opening(m1, footprint=se)
        m2 = morphology.binary_opening(m2, footprint=se)
        m3 = morphology.binary_opening(m3, footprint=se)

    return m1.astype(bool), m2.astype(bool), m3.astype(bool)


# ===================== Helpers: measurement =====================
def measure_gamma_pp(mask, px_nm_x, px_nm_y):
    """Plate-like γʺ: minor/major axis + aspect ratio."""
    lab = measure.label(mask, connectivity=2)
    if lab.max() == 0:
        return pd.DataFrame(columns=["minor_nm", "major_nm", "aspect_ratio"])

    props = measure.regionprops_table(
        lab, properties=("minor_axis_length", "major_axis_length")
    )
    df = pd.DataFrame(props).rename(columns={
        "minor_axis_length": "minor_px",
        "major_axis_length": "major_px"
    })

    scale_nm = np.sqrt(px_nm_x * px_nm_y)
    df["minor_nm"] = df["minor_px"] * scale_nm
    df["major_nm"] = df["major_px"] * scale_nm
    df["aspect_ratio"] = df["major_px"] / np.clip(df["minor_px"], 1e-12, None)

    return df[["minor_nm", "major_nm", "aspect_ratio"]]


def measure_gamma_p(mask, px_nm_x, px_nm_y):
    """Roughly spherical γ': equivalent radius from area."""
    lab = measure.label(mask, connectivity=2)
    if lab.max() == 0:
        return pd.DataFrame(columns=["radius_nm"])

    props = measure.regionprops_table(lab, properties=("area",))
    df = pd.DataFrame(props).rename(columns={"area": "area_px"})

    area_nm2_per_px = px_nm_x * px_nm_y
    df["radius_nm"] = np.sqrt((df["area_px"] * area_nm2_per_px) / np.pi)

    return df[["radius_nm"]]


# ===================== File discovery =====================
def find_series(root_dir: Path, master_glob: str):
    """
    Find main Exodus file and (optionally) companion -s### files.

    In your current setup: only one file like CS_2_out.e with internal time steps.
    """
    masters = [
        p for p in root_dir.glob(master_glob)
        if re.search(r"-s\d+$", p.name) is None
    ]
    companions = [
        p for p in root_dir.iterdir()
        if p.is_file() and re.search(r"\.(?:e|exo)-s\d+$", p.name, flags=re.I)
    ]
    companions.sort(key=lambda p: parse_step_from_name(p.name))
    return masters, companions


def read_companion(path: Path):
    """Read a companion .e-s### file (legacy)."""
    for forced in (".e", ".exo"):
        try:
            return pv.read(str(path), force_ext=forced)
        except Exception:
            pass
    return pv.read(str(path))


# ===================== Summaries =====================
def summarize_pp(per_step_rows):
    cols = [
        "step_index", "time_s", "n",
        "mean_minor_nm", "median_minor_nm",
        "mean_major_nm", "median_major_nm",
        "mean_AR", "median_AR"
    ]
    out = []
    for s, rec in sorted(per_step_rows.items()):
        mins = np.array([r[0] for r in rec], float)
        maxs = np.array([r[1] for r in rec], float)
        ars = np.array([r[2] for r in rec], float)
        times = [r[3] for r in rec if r[3] is not None]
        tmed = float(np.median(times)) if times else np.nan
        out.append([
            s, tmed, len(mins),
            mins.mean() if len(mins) else np.nan,
            np.median(mins) if len(mins) else np.nan,
            maxs.mean() if len(maxs) else np.nan,
            np.median(maxs) if len(maxs) else np.nan,
            ars.mean() if len(ars) else np.nan,
            np.median(ars) if len(ars) else np.nan,
        ])
    return pd.DataFrame(out, columns=cols)


def summarize_p(per_step_rows):
    cols = ["step_index", "time_s", "n", "mean_radius_nm", "median_radius_nm"]
    out = []
    for s, rec in sorted(per_step_rows.items()):
        rads = np.array([r[0] for r in rec], float)
        times = [r[1] for r in rec if r[1] is not None]
        tmed = float(np.median(times)) if times else np.nan
        out.append([
            s, tmed, len(rads),
            rads.mean() if len(rads) else np.nan,
            np.median(rads) if len(rads) else np.nan,
        ])
    return pd.DataFrame(out, columns=cols)


# ===================== Main processing =====================
masters, companions = find_series(root, args.master_glob)
if not masters and not companions:
    raise FileNotFoundError(f"No Exodus files found in {root}")

main_e = masters[0] if masters else None
reader = pv.get_reader(str(main_e)) if main_e else None
times = getattr(reader, "time_values", []) if reader is not None else []

# --------- Decide output path based on Exodus file name ---------
if args.out:
    out_path = args.out
else:
    if main_e is not None:
        # e.g. CS_2_out.e  ->  CS_2_out_AR.xlsx (same directory)
        stem = main_e.stem
        out_path = main_e.with_name(f"{stem}_AR.xlsx")
    else:
        # Fallback (very unlikely to be used)
        out_path = root / "AR_results.xlsx"
# ---------------------------------------------------------------

# Filter companions (mostly legacy, probably none for you now)
if args.file_stride > 1 and companions:
    companions = companions[::args.file_stride]
if args.max_files is not None and companions:
    companions = companions[:args.max_files]

# Determine sampling grid
if args.auto_res:
    if reader is not None and len(times) > 0:
        if hasattr(reader, "set_active_time_index"):
            reader.set_active_time_index(0)
        ds0 = reader.read()
    else:
        ds0 = read_companion(companions[0] if companions else masters[0])
    probe = ds0[0] if isinstance(ds0, pv.MultiBlock) else ds0
    NX, NY = choose_resolution_from_bounds(
        probe.bounds, args.unit_to_nm, args.L_MIN_NM
    )
else:
    NX, NY = args.nx, args.ny

print(f"Sampling grid: NX={NX}, NY={NY}")

# accumulators: per-step rows for ALL and PASSED
pp_all = {}     # step -> list of (minor_nm, major_nm, AR, time)
pp_pass = {}
p_all = {}      # step -> list of (radius_nm, time)
p_pass = {}


def add_pp(step, time_s, df_pp, passed_floor_nm):
    if len(df_pp):
        pp_all.setdefault(step, []).extend([
            (r.minor_nm, r.major_nm, r.aspect_ratio, time_s)
            for r in df_pp.itertuples()
        ])
        sel = df_pp[df_pp["minor_nm"] >= passed_floor_nm]
        if len(sel):
            pp_pass.setdefault(step, []).extend([
                (r.minor_nm, r.major_nm, r.aspect_ratio, time_s)
                for r in sel.itertuples()
            ])


def add_p(step, time_s, df_p, radius_floor_nm):
    if len(df_p):
        p_all.setdefault(step, []).extend([
            (r.radius_nm, time_s) for r in df_p.itertuples()
        ])
        sel = df_p[df_p["radius_nm"] >= radius_floor_nm]
        if len(sel):
            p_pass.setdefault(step, []).extend([
                (r.radius_nm, time_s) for r in sel.itertuples()
            ])


def process_dataset(ds, step_index, time_value):
    arrays = detect_eta_arrays_any(ds)
    info1 = arrays.get("eta_pv1")  # γʺ v1 (eta_gpp1)
    info2 = arrays.get("eta_pv2")  # γ'    (eta_gp)
    info3 = arrays.get("eta_pv3")  # γʺ v2 (eta_gpp2)

    if not any([info1, info2, info3]):
        return

    # choose a reference block for bounds (prefer γ', else γʺ)
    ref = info2 or info1 or info3
    ref_blk = get_block(ds, ref["block"])
    ref_arr, ref_loc = ref["name"], ref["location"]

    # sample reference once to set grid & extract bounds
    img_ref, dx, dy, _ = resample_to_uniform(
        ref_blk, ref_arr, ref_loc, NX, NY, args.unit_to_nm
    )

    # sample each variant on SAME grid
    def sample_var(info):
        if not info:
            return None
        blk = get_block(ds, info["block"])
        return resample_to_uniform(
            blk, info["name"], info["location"], NX, NY, args.unit_to_nm
        )[0]

    img1 = sample_var(info1)
    img2 = sample_var(info2)
    img3 = sample_var(info3)

    # If a variant missing, replace with -inf so it never wins argmax
    if img1 is None:
        img1 = np.full_like(img_ref, -np.inf, dtype=float)
    if img2 is None:
        img2 = np.full_like(img_ref, -np.inf, dtype=float)
    if img3 is None:
        img3 = np.full_like(img_ref, -np.inf, dtype=float)

    # exclusive segmentation (γʺ1, γ', γʺ2)
    m1, m2, m3 = exclusive_phase_masks(
        img1, img2, img3, dx, dy,
        tau1=args.tau_gpp1,
        tau2=args.tau_gp,
        tau3=args.tau_gpp2,
        L_MIN_NM=args.L_MIN_NM,
        open_frac=args.open_frac,
    )

    # measure
    df1 = measure_gamma_pp(m1, dx, dy)
    df3 = measure_gamma_pp(m3, dx, dy)
    df2 = measure_gamma_p(m2, dx, dy)

    # accumulate per-step rows (γʺ combines v1 + v2)
    if step_index >= args.start_step:
        tval = time_value
        if len(df1):
            add_pp(step_index, tval, df1, args.L_MIN_NM)
        if len(df3):
            add_pp(step_index, tval, df3, args.L_MIN_NM)
        if len(df2):
            add_p(step_index, tval, df2, args.R_MIN_NM)


# ---- process master internal steps (your main time series) ----
if reader is not None and len(times) > 0:
    for i, t in enumerate(times):
        if hasattr(reader, "set_active_time_value"):
            reader.set_active_time_value(float(t))
        else:
            reader.set_active_time_index(int(i))
        ds = reader.read()
        tval = read_time_value_any(ds)
        process_dataset(ds, i, tval if tval is not None else t)

# ---- process companion files (legacy) ----
for f in companions:
    ds = read_companion(f)
    step = parse_step_from_name(f.name)
    tval = read_time_value_any(ds)
    process_dataset(ds, step, tval)

# ===================== Build summaries & save =====================
gamma_pp_ALL = summarize_pp(pp_all)
gamma_pp_PASSED = summarize_pp(pp_pass)
gamma_p_ALL = summarize_p(p_all)
gamma_p_PASSED = summarize_p(p_pass)

print(
    f"Steps kept (>= {args.start_step}): "
    f"γʺ ALL={len(gamma_pp_ALL)}, γʺ PASSED={len(gamma_pp_PASSED)}, "
    f"γ' ALL={len(gamma_p_ALL)}, γ' PASSED={len(gamma_p_PASSED)}"
)

# Try Excel; if openpyxl missing, write CSVs
try:
    import openpyxl  # noqa: F401
    with pd.ExcelWriter(out_path) as xw:
        gamma_pp_ALL.to_excel(xw, sheet_name="gamma_pp_ALL", index=False)
        gamma_pp_PASSED.to_excel(xw, sheet_name="gamma_pp_PASSED", index=False)
        gamma_p_ALL.to_excel(xw, sheet_name="gamma_p_ALL", index=False)
        gamma_p_PASSED.to_excel(xw, sheet_name="gamma_p_PASSED", index=False)
    print("Saved:", out_path)
except Exception:
    print("openpyxl not available or Excel write failed; writing CSVs next to root instead.")
    base = out_path.with_suffix("")
    gamma_pp_ALL.to_csv(base.with_name(base.name + "_gamma_pp_ALL.csv"), index=False)
    gamma_pp_PASSED.to_csv(base.with_name(base.name + "_gamma_pp_PASSED.csv"), index=False)
    gamma_p_ALL.to_csv(base.with_name(base.name + "_gamma_p_ALL.csv"), index=False)
    gamma_p_PASSED.to_csv(base.with_name(base.name + "_gamma_p_PASSED.csv"), index=False)
    print("Saved CSVs with prefix:", base)
