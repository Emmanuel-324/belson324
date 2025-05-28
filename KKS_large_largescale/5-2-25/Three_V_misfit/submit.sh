#!/bin/bash                                                                     
#SBATCH --job-name   Three_V_misfit.i                                            
#SBATCH --out  Three_V_misfit.out        
#SBATCH --nodes 1                                                              
#SBATCH --ntasks 8
#SBATCH --cpus-per-task 4
#SBATCH --account amcorrosion                                                  
##SBATCH --partition normal_q                                                   
#SBATCH --time=120:00:00    
##SBATCH --mem=gpu:pascal:4                                                     
#SBATCH --export=NONE # this makes sure the compute environment is clean        
#SBATCH --mail-user emmanuel324@vt.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end


module reset
module load Miniforge3/24.1.2-0
source activate /home/emmanuel324/mambaforge3/envs/moose 
/home/emmanuel324/projects/belson324/belson324-opt  -i  Three_V_misfit.i







