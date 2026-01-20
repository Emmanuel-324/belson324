#!/bin/bash
#SBATCH --job-name=INB_2
#SBATCH --output=INB_2.%j.out
#SBATCH --error=INB_2.%j.err

#SBATCH --account=amcorrosion
#SBATCH --partition=normal_q
#SBATCH --time=5-00:00:00

#SBATCH --nodes=4
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=2

#SBATCH --mail-user=emmanuel324@vt.edu
#SBATCH --mail-type=BEGIN,END

# Clean module environment
module reset
module load EasyBuild/5.0.0
module load Miniforge3/24.11.3-0

# Activate your MOOSE env
source activate /home/emmanuel324/mambaforge3/envs/moose
export PATH="$CONDA_PREFIX/bin:$PATH"

# Go to the submission directory
cd "${SLURM_SUBMIT_DIR}"

# Make sure threading matches cpus-per-task
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# Run (do NOT hard-code -n; Slurm derives it from your SBATCH requests)
srun --export=ALL --mpi=pmix /home/emmanuel324/projects/belson324/belson324-opt -i INB_2.i
