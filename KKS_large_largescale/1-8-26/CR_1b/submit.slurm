#!/bin/bash
#SBATCH -A amcorrosion
#SBATCH -p normal_q
#SBATCH -t 06:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -J aspect_es01
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

# Use the Python from the env that has pyvista, scikit-image, pandas, etc.
# If unsure, run `which python` inside that env and paste the path here.
PY="/home/emmanuel324/mambaforge3/bin/python"

# Threading hints
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export VTK_SMP_MAX_THREADS=${SLURM_CPUS_PER_TASK}

"$PY" AR.py \
  --root "/home/emmanuel324/projects/belson324/KKS_large_largescale/1-8-26/CR_1b" \
  --master_glob "CR_1b_out.e" \
  --start_step 70 \
  --auto_res \
  --tau_gp 0.10 \
  --tau_gpp1 0.10 \
  --tau_gpp2 0.10 \
  --open_frac 0.25 \
  --L_MIN_NM 3.0 \
  --R_MIN_NM 2.0


