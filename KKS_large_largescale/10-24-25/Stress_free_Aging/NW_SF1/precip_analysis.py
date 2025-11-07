#!/usr/bin/env python
# --------------------------------------------------------------
# precip_analysis_nojoblib.py
# --------------------------------------------------------------
import argparse
from pathlib import Path
import numpy as np
from netCDF4 import Dataset
from scipy.ndimage import label, find_objects, center_of_mass
import pandas as pd
import matplotlib.pyplot as plt

# --------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Precipitate area & Hu invariants (no joblib)")
    p.add_argument("exodus_file", type=Path, help="Path to .e/.exo file")
    p.add_argument("--variants", nargs="+", default=["eta_pv1","eta_pv2","eta_pv3"],
                   help="Order-parameter names")
    p.add_argument("--threshold", type=float, default=0.5)
    p.add_argument("--outdir", type=Path, default=Path("precip_results"))
    p.add_argument("--plot-every", type=int, default=0,
                   help="Plot every N steps (0 = no plots)")
    return p.parse_args()

# --------------------------------------------------------------
def element_area(x, y):
    dx = np.min(np.diff(np.unique(x)))
    dy = np.min(np.diff(np.unique(y)))
    if not np.isclose(dx, dy):
        raise ValueError("Grid must be square and uniform")
    return 4 * dx * dx                     # quad-4 element area

def compute_moments(binary, centroid):
    yy, xx = np.indices(binary.shape)
    xx = xx - centroid[1]
    yy = yy - centroid[0]
    area = np.sum(binary)
    mu20 = np.sum(xx**2 * binary)
    mu02 = np.sum(yy**2 * binary)
    mu11 = np.sum(xx * yy * binary)
    return area, mu20, mu02, mu11

def hu_invariants(mu20, mu02, mu11):
    I1 = mu20 + mu02
    I2 = (mu20 - mu02)**2 + 4*mu11**2
    return I1, I2

# --------------------------------------------------------------
def process_step(step, nc, var, conn, elem_area, thresh, nx, ny):
    eta_nodes = nc.variables[var][step, :]
    eta_elem  = np.mean(eta_nodes[conn], axis=1)
    precip    = (eta_elem > thresh).astype(int)

    total_area = np.sum(precip) * elem_area

    mask_grid = precip.reshape(ny, nx)
    labeled, n_labels = label(mask_grid)
    objs = find_objects(labeled)

    invs = []
    for i, sl in enumerate(objs, 1):
        sub = mask_grid[sl]
        com = center_of_mass(sub)
        area, mu20, mu02, mu11 = compute_moments(sub, com)
        I1, I2 = hu_invariants(mu20, mu02, mu11)
        invs.append((i, area, I1, I2))

    return total_area, n_labels, invs

# --------------------------------------------------------------
def main():
    global args
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    with Dataset(args.exodus_file, 'r') as nc:
        x = nc.variables['coordx'][:]
        y = nc.variables['coordy'][:]
        conn = nc.variables['elem_conn'][:] - 1
        n_elem = conn.shape[0]
        elem_area = element_area(x, y)

        nx = int(np.sqrt(n_elem))
        ny = nx
        if nx*ny != n_elem:
            raise ValueError("Element count is not a perfect square – non-rectangular map?")

        n_steps = nc.dimensions['time_step'].size
        times = nc.variables['time_whole'][:] if 'time_whole' in nc.variables else np.arange(n_steps)

        results = {v: {"time": [], "total_area": [], "n_part": [], "particles": []} for v in args.variants}

        for step in range(n_steps):
            if step % 20 == 0:
                print(f"Step {step}/{n_steps-1}")

            for var in args.variants:
                total_area, n_part, invs = process_step(
                    step, nc, var, conn, elem_area, args.threshold, nx, ny
                )
                results[var]["time"].append(times[step])
                results[var]["total_area"].append(total_area)
                results[var]["n_part"].append(n_part)
                results[var]["particles"].append(invs)

            # ---- optional plot ----
            if args.plot_every and (step % args.plot_every == 0):
                plt.figure(figsize=(4*len(args.variants), 4))
                for i, var in enumerate(args.variants, 1):
                    grid = nc.variables[var][step, :].reshape(ny, nx)
                    plt.subplot(1, len(args.variants), i)
                    im = plt.imshow(grid, cmap='viridis')
                    plt.title(f"{var} – step {step}")
                    plt.colorbar(im, fraction=0.046, pad=0.04)
                plt.tight_layout()
                plt.savefig(args.outdir / f"step_{step:05d}.png")
                plt.close()

        # ---- write CSV ----
        for var in args.variants:
            df = pd.DataFrame({
                "time": results[var]["time"],
                "total_area_m2": results[var]["total_area"],
                "n_particles": results[var]["n_part"]
            })
            df.to_csv(args.outdir / f"{var}_summary.csv", index=False)

            rows = []
            for t, invs in zip(results[var]["time"], results[var]["particles"]):
                for pid, area, I1, I2 in invs:
                    rows.append([t, pid, area, I1, I2])
            pd.DataFrame(rows, columns=["time","particle_id","area_m2","Hu_I1","Hu_I2"])\
              .to_csv(args.outdir / f"{var}_particles.csv", index=False)

    print(f"\nAll done – results in: {args.outdir.resolve()}")

if __name__ == "__main__":
    main()