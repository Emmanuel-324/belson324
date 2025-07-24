#!/bin/bash                                                                     
#SBATCH --job-name=grain_growth_2D_random.i                                            
#SBATCH --output=grain_growth_2D_random.%j.out 
#SBATCH --error=grain_growth_2D_random.%j.err       
#SBATCH --nodes=1                                                              
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=8
#SBATCH --account=amcorrosion                                                  
##SBATCH --partition=normal_q                                                   
#SBATCH --time=120:00:00    
#SBATCH --export=NONE # this makes sure the compute environment is clean        
#SBATCH --mail-user emmanuel324@vt.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end


module reset
module load EasyBuild/5.0.0
module load Miniforge3/24.11.3-0
source activate /home/emmanuel324/mambaforge3/envs/moose 

export PATH="$CONDA_PREFIX/bin:$PATH"

# Go to job submission directory
cd "${SLURM_SUBMIT_DIR}"

srun --export=ALL --mpi=pmi2 -n 4 /home/emmanuel324/projects/belson324/belson324-opt  -i  grain_growth_2D_random.i







