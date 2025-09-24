# ====== Aspect ratio from Exodus (.e) per timestep & variant ======
from pathlib import Path
import re
import numpy as np
import pandas as pd
import pyvista as pv
from skimage import filters, measure, morphology

# ---------------- CONFIG (edit these) ----------------
root = Path(r"/home/emmanuel324/projects/belson324/KKS_large_largescale/9-12-25/ES_02")  # folder containing .e and .e-s###
MASTER_GLOB   = "*.e"                       # master .e pattern (e.g., 'Sing1_adt_out.e')
VAR_NAMES     = ["eta_pv1", "eta_pv2", "eta_pv3"]  # what to look for (with fallbacks)
NX, NY        = 1024, 1024                  # raster resolution (uniform XY grid)
THRESHOLD     = 0.5                         # or "otsu" for adaptive
MIN_AREA_PX   = 20                          # ignore specks smaller than this (in pixels)
FILE_STRIDE   = 1                           # process every k-th .e-s### file (1 = all)
MAX_FILES     = None                        # limit count when testing (e.g., 50); None = all
UNIT_TO_NM    = 1.0                         # multiply Exodus coordinates by this to convert to nm (1.0 if already nm)

# ---- PyVista compatibility: UniformGrid vs ImageData ----
UniformGrid = getattr(pv, "UniformGrid", pv.ImageData)

# ---------------- Helpers ----------------
def parse_step_from_name(name: str) -> int:
    """'file.e' -> 0, 'file.e-s002' -> 2, 'file.e-s505' -> 505"""
    m = re.search(r"\.e(?:-s(\d+))?$", name)
    return int(m.group(1)) if (m and m.group(1)) else 0

def read_time_value_any(ds) -> float | None:
    """Try to read a physical time scalar from field_data at root or leaves."""
    def grab(fd):
        for key in ("TimeValue", "time", "TIME", "Time"):
            if key in fd:
                v = fd[key]
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
    """
    Return (actual_name, location) where location is 'point' or 'cell', or (None, None) if not found.
    Searches point_data first (nodal), then cell_data.
    accepted_names: list of lowercase names to try, e.g., ['eta_pv1','eta1']
    """
    pkeys = list(dset.point_data.keys())
    ckeys = list(dset.cell_data.keys())
    lower_p = {k.lower(): k for k in pkeys}
    lower_c = {k.lower(): k for k in ckeys}

    for candidate in accepted_names:
        if candidate in lower_p:
            return lower_p[candidate], "point"
    for candidate in accepted_names:
        if candidate in lower_c:
            return lower_c[candidate], "cell"
    return None, None

def detect_eta_arrays_any(ds):
    """
    Works for single datasets or MultiBlock.
    Returns dict variant-> {'name': str, 'location': 'point'|'cell', 'block': block_index_tuple}
    Picks the first block where each eta is found.
    """
    want = {
        "eta_pv1": ["eta_pv1", "eta1"],
        "eta_pv2": ["eta_pv2", "eta2"],
        "eta_pv3": ["eta_pv3", "eta3"],
    }
    results = {k: None for k in want}

    def maybe_set(tag, val, block_idx):
        if val is None or results[tag] is not None:
            return
        name, loc = val
        if name is None:
            return
        results[tag] = {"name": name, "location": loc, "block": block_idx}

    if isinstance(ds, pv.MultiBlock):
        for bi, child in enumerate(ds):
            if child is None:
                continue
            if isinstance(child, pv.MultiBlock):
                for bj, gchild in enumerate(child):
                    if gchild is None:
                        continue
                    for tag, alts in want.items():
                        val = find_variant_in_dataset(gchild, alts)
                        maybe_set(tag, val, (bi, bj))
            else:
                for tag, alts in want.items():
                    val = find_variant_in_dataset(child, alts)
                    maybe_set(tag, val, (bi,))
    else:
        for tag, alts in want.items():
            name, loc = find_variant_in_dataset(ds, alts)
            if name:
                results[tag] = {"name": name, "location": loc, "block": ()}
    return results

def get_block(ds, block_idx):
    """Return the block by tuple index; () means ds itself."""
    if not block_idx:
        return ds
    out = ds
    for i in block_idx:
        out = out[i]
    return out

def resample_to_uniform(block, array_name, location, nx=1024, ny=1024):
    """
    Sample the dataset onto a uniform XY grid. Returns (img, px_nm_x, px_nm_y).
    Assumes quasi-2D (thin z); uses mid-plane z.
    """
    # Convert cell -> point if needed
    if location == "cell":
        block = block.cell_data_to_point_data()

    # If volume, extract surface helps; if already PolyData, it's cheap
    try:
        surf = block.extract_surface()
        work = surf if surf.n_points > 0 else block
    except Exception:
        work = block

    xmin, xmax, ymin, ymax, zmin, zmax = work.bounds
    # convert units to nm if needed
    xmin *= UNIT_TO_NM; xmax *= UNIT_TO_NM
    ymin *= UNIT_TO_NM; ymax *= UNIT_TO_NM
    zmin *= UNIT_TO_NM; zmax *= UNIT_TO_NM

    if nx < 2 or ny < 2:
        raise ValueError("nx, ny must be >= 2")
    z0 = 0.5 * (zmin + zmax)
    dx = (xmax - xmin) / (nx - 1)
    dy = (ymax - ymin) / (ny - 1)

    grid = UniformGrid()
    grid.dimensions = (nx, ny, 1)
    grid.origin     = (xmin, ymin, z0)
    grid.spacing    = (dx,   dy,   max(zmax - zmin, 1.0))

    sampled = grid.sample(work)
    if array_name not in sampled.point_data:
        raise KeyError(f"Array '{array_name}' not found after sampling.")
    vals = np.asarray(sampled.point_data[array_name])
    if vals.size != nx * ny:
        raise RuntimeError(f"Sampled array size mismatch: {vals.size} != {nx*ny}")
    img = vals.reshape((ny, nx))
    px_nm_x = dx
    px_nm_y = dy
    return img, px_nm_x, px_nm_y

def mask_from_eta_image(img, thresh=0.5, min_area_px=0):
    arr = np.asarray(img, float)
    if isinstance(thresh, str) and thresh.lower() == "otsu":
        t = filters.threshold_otsu(arr)
    else:
        t = float(thresh)
    mask = arr >= t
    if min_area_px and min_area_px > 0:
        mask = morphology.remove_small_objects(mask.astype(bool), min_size=min_area_px)
    return mask

def measure_aspect_ratios(mask, px_nm_x=1.0, px_nm_y=1.0):
    labeled = measure.label(mask, connectivity=2)
    if labeled.max() == 0:
        return pd.DataFrame(columns=[
            "label","area_px","area_nm2","centroid_row","centroid_col",
            "major_px","minor_px","major_nm","minor_nm","orientation_deg","aspect_ratio"
        ])
    props = measure.regionprops_table(
        labeled,
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
    # Convert to nm (geometric mean if anisotropic pixels; exact if square pixels)
    px_to_nm = np.sqrt(px_nm_x * px_nm_y)
    df["major_nm"] = df["major_px"] * px_to_nm
    df["minor_nm"] = df["minor_px"] * px_to_nm
    df["area_nm2"] = df["area_px"] * (px_nm_x * px_nm_y)
    df["orientation_deg"] = np.degrees(df["orientation_rad"])
    eps = 1e-12
    df["aspect_ratio"] = df["major_px"] / np.maximum(df["minor_px"], eps)
    return df[["label","area_px","area_nm2","centroid_row","centroid_col",
               "major_px","minor_px","major_nm","minor_nm","orientation_deg","aspect_ratio"]]

# ---------------- Collect files ----------------
masters = [p for p in root.glob(MASTER_GLOB) if "-s" not in p.name]
if not masters:
    raise FileNotFoundError("No master .e file found.")
main_e = masters[0]
reader = pv.get_reader(str(main_e))
times = getattr(reader, "time_values", [])

companions = sorted(root.glob(main_e.name + "-s*"), key=lambda p: parse_step_from_name(p.name))
if FILE_STRIDE > 1:
    companions = companions[::FILE_STRIDE]
if MAX_FILES is not None:
    companions = companions[:MAX_FILES]

# ---------------- Process ----------------
rows = []

# 1) master internal steps (if any)
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
        block = get_block(ds, info["block"])
        arr_name, loc = info["name"], info["location"]

        img, px_x, px_y = resample_to_uniform(block, arr_name, loc, NX, NY)
        mask = mask_from_eta_image(img, THRESHOLD, MIN_AREA_PX)
        dfm = measure_aspect_ratios(mask, px_x, px_y)
        dfm.insert(0, "variant", var)
        dfm.insert(0, "time_value", tval if tval is not None else t)
        dfm.insert(0, "step_index", i)
        dfm.insert(0, "file", main_e.name)
        rows.append(dfm)

# 2) companion .e-s### files (bulk of steps)
for f in companions:
    ds = pv.read(str(f), force_ext=".e")  # force Exodus reader for .e-s### files
    step = parse_step_from_name(f.name)
    tval = read_time_value_any(ds)

    arrays = detect_eta_arrays_any(ds)
    for var in VAR_NAMES:
        info = arrays.get(var)
        if not info:
            continue
        block = get_block(ds, info["block"])
        arr_name, loc = info["name"], info["location"]

        img, px_x, px_y = resample_to_uniform(block, arr_name, loc, NX, NY)
        mask = mask_from_eta_image(img, THRESHOLD, MIN_AREA_PX)
        dfm = measure_aspect_ratios(mask, px_x, px_y)
        dfm.insert(0, "variant", var)
        dfm.insert(0, "time_value", tval if tval is not None else np.nan)
        dfm.insert(0, "step_index", step)
        dfm.insert(0, "file", f.name)
        rows.append(dfm)

# ---------------- Save ----------------
out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
    columns=["file","step_index","time_value","variant","label","area_px","area_nm2",
             "centroid_row","centroid_col","major_px","minor_px","major_nm","minor_nm",
             "orientation_deg","aspect_ratio"]
)
out_path = root / "aspect_ratios_by_precip.csv"
out.to_csv(out_path, index=False)
print(f"Rows: {len(out)}")
print("Saved:", out_path)
display(out.head(10))


