#!/bin/bash                                                                     
#SBATCH --job-name=BiC_2a.i                                            
#SBATCH --output=BiC_2a.%j.out 
#SBATCH --error=BiC_2a.%j.err       
#SBATCH --nodes=1                                                              
#SBATCH --ntasks-per-node=34
#SBATCH --cpus-per-task=2
#SBATCH --account=amcorrosion                                                  
#SBATCH --partition=normal_q                                                   
#SBATCH --time=7-00:00:00     
##SBATCH --mem=gpu:pascal:4
##SBATCH --mem=254G                                                     
##SBATCH --export=NONE # this makes sure the compute environment is clean        
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

srun --export=ALL --mpi=pmi2 -n 34 /home/emmanuel324/projects/belson324/belson324-opt  -i  BiC_2a.i





