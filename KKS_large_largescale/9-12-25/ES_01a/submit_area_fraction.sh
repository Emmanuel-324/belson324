#!/bin/bash
#SBATCH -A amcorrosion
#SBATCH -p normal_q
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH -J area_frac
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

# --- EDIT THESE ARGS (define early to avoid unbound variable errors) ---
EXO="/home/emmanuel324/projects/belson324/KKS_large_largescale/9-12-25/ES_01a/ES_01a_out.e"
VARS="h_pv1_aux h_pv2_aux h_pv3_aux h_m_aux"  # Updated to match MOOSE AuxVariables
OUT=""  # optional: set to a full path for CSV; leave empty to auto-name next to the .e

# Debugging output (now after definitions)
echo "Current working directory: $PWD"
ls -l
which python
echo "Python script path: /home/emmanuel324/projects/belson324/KKS_large_largescale/9-12-25/ES_01a/area_fraction_percent.py"
echo "Exodus file path: $EXO"
echo "Variables: $VARS"

# Use the Python from the env that has netCDF4 + pandas/numpy
PY="/home/emmanuel324/mambaforge3/bin/python"

# Threading hints
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Verify files exist
echo "Checking for Python script and .e file..."
ls -l /home/emmanuel324/projects/belson324/KKS_large_largescale/9-12-25/ES_01a/area_fraction_percent.py
ls -l "$EXO"

# Run
$PY /home/emmanuel324/projects/belson324/KKS_large_largescale/9-12-25/ES_01a/area_fraction_percent.py \
  --exo "$EXO" \
  --vars $VARS \
  ${OUT:+--out "$OUT"}