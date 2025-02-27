#!/bin/bash                                                                     
#SBATCH --job-name    slip_mode1_eigen                                      
#SBATCH --out  slip_mode1_eigen.out        
#SBATCH --nodes 1                                                              
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 8
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
module load OpenMPI
source activate /home/emmanuel324/mambaforge3/envs/moose 
mpirun -np 2 /home/emmanuel324/projects/belson324/belson324-opt  -i  slip_mode1_eigen.i







