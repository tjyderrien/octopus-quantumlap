#!/bin/bash  
#SBATCH --account=your-account-acc1
#SBATCH --partition=queue12h
#SBATCH -N 4 # number of nodes 
#SBATCH --clusters=eclipse
#SBATCH -t 0-2:00 # time (D-HH:MM) 
#SBATCH -o slurm.%N.%j.out # STDOUT 
#SBATCH -e slurm.%N.%j.err # STDERR 
#SBATCH --mail-type=END,FAIL # notifications for job done & fail 
#SBATCH --mail-user=youremail@eli-beams.eu # send-to address

module load OpenMPI/1.8.4-GCC-4.8.4
mpirun -npernode 2  ./stream_mpi_test
echo "Done"


