#!/bin/bash                                                                     
#SBATCH --job-name update_method_test.i                                               
#SBATCH --out  update_method_test.out        
#SBATCH --nodes 1                                                               
#SBATCH --ntasks 2
#SBATCH --cpus-per-task 8
#SBATCH --account amcorrosion                                                  
##SBATCH --partition normal_q                                                   
#SBATCH --time=100:00:00    
##SBATCH --mem=gpu:pascal:4                                                     
#SBATCH --export=NONE # this makes sure the compute environment is clean        
#SBATCH --mail-user emmanuel324@vt.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end


module reset
module load apps site/tinkercliffs/easybuild/setup
module load apps site/tinkercliffs-rome/easybuild/setup
module load OpenMPI CMake
export CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77
mamba activate python-stack


cd "${submitdir}"
#rm kks_multiphase_largescale_large_csv*
#rm kks_multiphase_largescale_large_out.e*
mpirun --mca orte_base_help_aggregate 0 -np 2 /home/emmanuel324/projects/belson/belson-opt  -i  update_method_test.i



