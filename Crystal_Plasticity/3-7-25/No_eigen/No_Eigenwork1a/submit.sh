#!/bin/bash                                                                     
#SBATCH --job-name   No_Eigenwork1a                                            
#SBATCH --out  No_Eigenwork1a.out        
#SBATCH --nodes 1                                                              
#SBATCH --ntasks-per-node 20
#SBATCH --cpus-per-task 4
#SBATCH --account amcorrosion                                                  
##SBATCH --partition normal_q                                                   
#SBATCH --time=130:00:00    
##SBATCH --mem=gpu:pascal:4                                                     
#SBATCH --export=NONE # this makes sure the compute environment is clean        
#SBATCH --mail-user emmanuel324@vt.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end


module reset
module load Miniforge3/24.1.2-0
source activate /home/emmanuel324/mambaforge3/envs/moose 
/home/emmanuel324/projects/belson324/belson324-opt  -i  No_Eigenwork1a.i







