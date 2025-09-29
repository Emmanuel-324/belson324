#!/bin/bash
#SBATCH -A amcorrosion
#SBATCH -p normal_q
#SBATCH -t 06:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -J aspect_single
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export VTK_SMP_MAX_THREADS=${SLURM_CPUS_PER_TASK}
python aspect_ratio_new.py
SLURM

chmod +x submit.slurm
sbatch submit.slurm
