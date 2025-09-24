from pathlib import Path
import re
import pyvista as pv

# --- Point to your folder (use raw string r'...' or forward slashes) ---
root = Path(r"/home/emmanuel324/projects/belson324/KKS_large_largescale/9-12-25/ES_01")

# Sanity check: does the directory exist?
assert root.exists(), f"Directory not found: {root}"
print("Directory exists:", root)
print("Sample files:", [p.name for p in list(root.iterdir())[:8]])

# ---- Helpers ----
def parse_step_from_name(name: str) -> int:
    m = re.search(r"\.(?:e|exo)-s(\d+)$", name, flags=re.I)
    return int(m.group(1)) if m else 0

def find_series(root: Path):
    # masters: *.e or *.exo that do NOT have -s### in the name
    masters = []
    for pat in ("*.e", "*.E", "*.exo", "*.EXO"):
        masters.extend([p for p in root.glob(pat) if re.search(r"-s\d+$", p.name, flags=re.I) is None])

    # companions: .e-s### or .exo-s###
    companions = [p for p in root.iterdir()
                  if p.is_file() and re.search(r"\.(?:e|exo)-s\d+$", p.name, flags=re.I)]
    companions.sort(key=lambda p: parse_step_from_name(p.name))
    return masters, companions

def read_companion(path: Path):
    # Force the Exodus reader for weird suffixes like .e-s004
    for forced in (".e", ".exo"):
        try:
            return pv.read(str(path), force_ext=forced)
        except Exception:
            pass
    return pv.read(str(path))

# ---- Collect
masters, companions = find_series(root)
print(f"Masters found: {[m.name for m in masters]}")
print(f"Companions found: {len(companions)} (first 5: {[p.name for p in companions[:5]]})")

# ---- Try reading master time-series if present
if masters:
    main_e = masters[0]
    reader = pv.get_reader(str(main_e))
    times = getattr(reader, "time_values", [])
    print(f"Using master: {main_e.name} | internal steps: {len(times)}")

    # Read a couple of internal steps
    for i, t in enumerate(times[:2]):
        if hasattr(reader, "set_active_time_value"):
            reader.set_active_time_value(float(t))
        else:
            reader.set_active_time_index(int(i))
        ds = reader.read()
        nblocks = ds.n_blocks if isinstance(ds, pv.MultiBlock) else 1
        print(f"[master] i={i} time={t} blocks={nblocks}")

# ---- Always test a few companion files too
for f in companions[:3]:
    ds = read_companion(f)
    nblocks = ds.n_blocks if isinstance(ds, pv.MultiBlock) else 1
    print(f"[companion] {f.name} step={parse_step_from_name(f.name)} blocks={nblocks}")
