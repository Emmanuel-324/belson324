from pathlib import Path
import re
import numpy as np
import pandas as pd
import pyvista as pv
from skimage import filters, measure, morphology

# ================== CONFIG ==================
root = Path(r"/home/emmanuel324/projects/belson324/KKS_large_largescale/9-12-25/EB_02")

VAR_NAMES = ["eta_pv1", "eta_pv2", "eta_pv3"]  # search order
THRESHOLD = "otsu"                              # or a float like 0.5
UNIT_TO_NM = 1.0                                # multiply coordinates by this to get nm
L_MIN_NM = 3.0                                  # smallest minor axis (nm) you want considered reliable
AUTO_RESOLVE = True                             # auto-pick NX/NY from L_MIN_NM and domain size
NX, NY = 1024, 1024                             # used if AUTO_RESOLVE=False
FILE_STRIDE = 1                                 # process every k-th companion file
MAX_FILES = None                                # e.g., 100 for a quick test; None = all
MASTER_GLOB = "*.e"                             # find master .e (no -s###)
OUT_CSV_NAME = "aspect_ratios_by_precip.csv"    # output file

# phase mapping
PHASE_MAP = {"eta_pv1": "gamma_pp", "eta_pv3": "gamma_pp", "eta_pv2": "gamma_p"}

# PyVista compat
UniformGrid = getattr(pv, "UniformGrid", pv.ImageData)

# ================== HELPERS ==================
def parse_step_from_name(name: str) -> int:
    m = re.search(r"\.(?:e|exo)-s(\d+)$", name, flags=re.I)
    return int(m.group(1)) if m else 0

def read_time_value_any(ds) -> float | None:
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
        for child in ds:
            if child is None: 
                continue
            if hasattr(child, "field_data"):
                t = grab(child.field_data)
                if t is not None:
                    return t
    return None

def find_variant_in_dataset(dset, accepted_names):
    pkeys = list(dset.point_data.keys())
    ckeys = list(dset.cell_data.keys())
    lp = {k.lower(): k for k in pkeys}
    lc = {k.lower(): k for k in ckeys}
    for nm in accepted_names:
        if nm in lp:
            return lp[nm], "point"
    for nm in accepted_names:
        if nm in lc:
            return lc[nm], "cell"
    return None, None

def detect_eta_arrays_any(ds):
    want = {
        "eta_pv1": ["eta_pv1", "eta1"],
        "eta_pv2": ["eta_pv2", "eta2"],
        "eta_pv3": ["eta_pv3", "eta3"],
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
    if not idx_tuple:
        return ds
    cur = ds
    for k in idx_tuple:
        cur = cur[k]
    return cur

def choose_resolution_from_bounds(bounds, unit_to_nm: float, l_min_nm: float, min_dim=256, max_dim=2048):
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    xmin *= unit_to_nm; xmax *= unit_to_nm
    ymin *= unit_to_nm; ymax *= unit_to_nm
    Lx = max(xmax - xmin, 1e-9)
    Ly = max(ymax - ymin, 1e-9)
    p_target = max(l_min_nm / 6.0, 1e-6)     # ~6 px across L_min
    NX = int(np.clip(round(Lx / p_target) + 1, min_dim, max_dim))
    NY = int(np.clip(round(Ly / p_target) + 1, min_dim, max_dim))
    return NX, NY

def resample_to_uniform(block, array_name, location, nx, ny, unit_to_nm):
    if location == "cell":
        block = block.cell_data_to_point_data()
    try:
        surf = block.extract_surface()
        work = surf if surf.n_points > 0 else block
    except Exception:
        work = block
    xmin, xmax, ymin, ymax, zmin, zmax = work.bounds
    xmin *= unit_to_nm; xmax *= unit_to_nm
    ymin *= unit_to_nm; ymax *= unit_to_nm
    zmin *= unit_to_nm; zmax *= unit_to_nm
    z0 = 0.5 * (zmin + zmax)
    dx = (xmax - xmin) / max(nx - 1, 1)
    dy = (ymax - ymin) / max(ny - 1, 1)
    grid = UniformGrid()
    grid.dimensions = (nx, ny, 1)
    grid.origin = (xmin, ymin, z0)
    grid.spacing = (dx, dy, max(zmax - zmin, 1.0))
    smp = grid.sample(work)
    if array_name not in smp.point_data:
        raise KeyError(f"Array '{array_name}' not found after sampling.")
    vals = np.asarray(smp.point_data[array_name])
    img = vals.reshape((ny, nx))
    return img, dx, dy

def make_mask(img, thresh, min_area_px: int):
    arr = np.asarray(img, float)
    if isinstance(thresh, str) and thresh.lower() == "otsu":
        t = filters.threshold_otsu(arr)
    else:
        t = float(thresh)
    m = arr >= t
    if min_area_px and min_area_px > 0:
        m = morphology.remove_small_objects(m.astype(bool), min_size=int(min_area_px))
    return m

def measure_gamma_pp(mask, px_nm_x, px_nm_y):
    lab = measure.label(mask, connectivity=2)
    if lab.max() == 0:
        return pd.DataFrame(columns=["label","area_px","area_nm2","centroid_row","centroid_col",
                                     "major_px","minor_px","major_nm","minor_nm","orientation_deg","aspect_ratio"])
    props = measure.regionprops_table(
        lab,
        properties=("label","area","centroid","orientation","major_axis_length","minor_axis_length")
    )
    df = pd.DataFrame(props).rename(columns={
        "area": "area_px",
        "centroid-0": "centroid_row",
        "centroid-1": "centroid_col",
        "orientation": "orientation_rad",
        "major_axis_length": "major_px",
        "minor_axis_length": "minor_px"
    })
    scale_nm = np.sqrt(px_nm_x * px_nm_y)
    df["area_nm2"] = df["area_px"] * (px_nm_x * px_nm_y)
    df["major_nm"] = df["major_px"] * scale_nm
    df["minor_nm"] = df["minor_px"] * scale_nm
    df["orientation_deg"] = np.degrees(df["orientation_rad"])
    df["aspect_ratio"] = df["major_px"] / np.clip(df["minor_px"], 1e-12, None)
    return df[["label","area_px","area_nm2","major_px","minor_px","major_nm","minor_nm","orientation_deg","aspect_ratio","centroid_row","centroid_col"]]

def measure_gamma_p_radius(mask, px_nm_x, px_nm_y):
    lab = measure.label(mask, connectivity=2)
    if lab.max() == 0:
        return pd.DataFrame(columns=["label","area_px","area_nm2","radius_nm","centroid_row","centroid_col"])
    props = measure.regionprops_table(
        lab,
        properties=("label","area","centroid")
    )
    df = pd.DataFrame(props).rename(columns={
        "area": "area_px",
        "centroid-0": "centroid_row",
        "centroid-1": "centroid_col",
    })
    area_nm2_per_px = (px_nm_x * px_nm_y)
    df["area_nm2"] = df["area_px"] * area_nm2_per_px
    df["radius_nm"] = np.sqrt(df["area_nm2"] / np.pi)
    return df[["label","area_px","area_nm2","radius_nm","centroid_row","centroid_col"]]

def find_series(root_dir: Path):
    masters = [p for p in root_dir.glob(MASTER_GLOB) if re.search(r"-s\d+$", p.name) is None]
    companions = [p for p in root_dir.iterdir()
                  if p.is_file() and re.search(r"\.(?:e|exo)-s\d+$", p.name, flags=re.I)]
    companions.sort(key=lambda p: parse_step_from_name(p.name))
    return masters, companions

def read_companion(path: Path):
    for forced in (".e", ".exo"):
        try:
            return pv.read(str(path), force_ext=forced)
        except Exception:
            pass
    return pv.read(str(path))

# ================== COLLECT ==================
masters, companions = find_series(root)
if not masters and not companions:
    raise FileNotFoundError(f"No Exodus files found in {root}")

main_e = masters[0] if masters else None
reader = pv.get_reader(str(main_e)) if main_e else None
times = getattr(reader, "time_values", []) if reader is not None else []

if FILE_STRIDE > 1 and companions:
    companions = companions[::FILE_STRIDE]
if MAX_FILES is not None and companions:
    companions = companions[:MAX_FILES]

# Determine NX/NY automatically (once) if requested
if AUTO_RESOLVE:
    # probe bounds from master or first companion
    if reader is not None:
        if hasattr(reader, "set_active_time_index"):
            reader.set_active_time_index(0)
        ds0 = reader.read()
    else:
        ds0 = read_companion(companions[0])
    probe = ds0[0] if isinstance(ds0, pv.MultiBlock) else ds0
    NX, NY = choose_resolution_from_bounds(probe.bounds, UNIT_TO_NM, L_MIN_NM)

rows = []

# ================== PROCESS: MASTER INTERNAL STEPS ==================
if reader is not None and len(times) > 0:
    for i, t in enumerate(times):
        if hasattr(reader, "set_active_time_value"):
            reader.set_active_time_value(float(t))
        else:
            reader.set_active_time_index(int(i))
        ds = reader.read()
        tval = read_time_value_any(ds)
        arrays = detect_eta_arrays_any(ds)
        for var in VAR_NAMES:
            info = arrays.get(var)
            if not info:
                continue
            blk = get_block(ds, info["block"])
            arr_name, loc = info["name"], info["location"]
            img, dx, dy = resample_to_uniform(blk, arr_name, loc, NX, NY, UNIT_TO_NM)
            px_geom = np.sqrt(dx * dy)
            if var in ("eta_pv1", "eta_pv3"):
                # dynamic min area from L_MIN_NM
                minor_px = max(L_MIN_NM / px_geom, 0.0)
                min_area_px = int(np.ceil((np.pi/4.0) * minor_px**2))
                m = make_mask(img, THRESHOLD, min_area_px)
                dfm = measure_gamma_pp(m, dx, dy)
                if len(dfm):
                    dfm.insert(0, "phase", PHASE_MAP[var])
                    dfm.insert(0, "variant", var)
                    dfm.insert(0, "time_value", tval if tval is not None else t)
                    dfm.insert(0, "step_index", i)
                    dfm.insert(0, "file", main_e.name)
                    dfm["passes_Lmin_nm"] = dfm["minor_nm"] >= L_MIN_NM
                    rows.append(dfm)
            else:  # eta_pv2 → gamma'
                m = make_mask(img, THRESHOLD, 0)
                dfm = measure_gamma_p_radius(m, dx, dy)
                if len(dfm):
                    dfm.insert(0, "phase", PHASE_MAP[var])
                    dfm.insert(0, "variant", var)
                    dfm.insert(0, "time_value", tval if tval is not None else t)
                    dfm.insert(0, "step_index", i)
                    dfm.insert(0, "file", main_e.name)
                    rows.append(dfm)

# ================== PROCESS: COMPANION FILES ==================
for f in companions:
    ds = read_companion(f)
    step = parse_step_from_name(f.name)
    tval = read_time_value_any(ds)
    arrays = detect_eta_arrays_any(ds)
    for var in VAR_NAMES:
        info = arrays.get(var)
        if not info:
            continue
        blk = get_block(ds, info["block"])
        arr_name, loc = info["name"], info["location"]
        img, dx, dy = resample_to_uniform(blk, arr_name, loc, NX, NY, UNIT_TO_NM)
        px_geom = np.sqrt(dx * dy)
        if var in ("eta_pv1", "eta_pv3"):
            minor_px = max(L_MIN_NM / px_geom, 0.0)
            min_area_px = int(np.ceil((np.pi/4.0) * minor_px**2))
            m = make_mask(img, THRESHOLD, min_area_px)
            dfm = measure_gamma_pp(m, dx, dy)
            if len(dfm):
                dfm.insert(0, "phase", PHASE_MAP[var])
                dfm.insert(0, "variant", var)
                dfm.insert(0, "time_value", tval if tval is not None else np.nan)
                dfm.insert(0, "step_index", step)
                dfm.insert(0, "file", f.name)
                dfm["passes_Lmin_nm"] = dfm["minor_nm"] >= L_MIN_NM
                rows.append(dfm)
        else:  # eta_pv2 → gamma'
            m = make_mask(img, THRESHOLD, 0)
            dfm = measure_gamma_p_radius(m, dx, dy)
            if len(dfm):
                dfm.insert(0, "phase", PHASE_MAP[var])
                dfm.insert(0, "variant", var)
                dfm.insert(0, "time_value", tval if tval is not None else np.nan)
                dfm.insert(0, "step_index", step)
                dfm.insert(0, "file", f.name)
                rows.append(dfm)

# ================== SAVE ==================
if rows:
    out = pd.concat(rows, ignore_index=True)
else:
    # unified schema with NaNs for missing columns
    out = pd.DataFrame(columns=[
        "file","step_index","time_value","variant","phase",
        "label","area_px","area_nm2",
        "major_px","minor_px","major_nm","minor_nm","orientation_deg","aspect_ratio",
        "radius_nm","centroid_row","centroid_col","passes_Lmin_nm"
    ])

# ensure all expected columns exist
for col in ["major_px","minor_px","major_nm","minor_nm","orientation_deg","aspect_ratio","radius_nm","passes_Lmin_nm"]:
    if col not in out.columns:
        out[col] = np.nan
out_path = root / OUT_CSV_NAME
out.to_csv(out_path, index=False)
print(f"NX, NY used: {NX}, {NY}")
print(f"Rows: {len(out)}")
print("Saved:", out_path)
