#######################################################
#example SGE batch script for a plain MPI program using
#512 cores (16 nodes x 32 cores per node) with 32 MPI 
#tasks per node
#executed under the intel runtime 
#######################################################

## run in /bin/bash
#$ -S /bin/bash
## do not join stdout and stderr
#$ -j n
## name of the job
#$ -N TutoBenzene
## execute job from the current working directory
#$ -cwd
## do not send mail
## $ -m n
## request 16 nodes (x 32 cores), must be a multiple of 32  
#$ -pe impi_hydra 32 
## run for 30 minutes
#$ -l h_rt=00:05:00

module load intel impi
module load octopus/5.0.0

## to save memory at high core counts, the "connectionless user datagram protocol" 
## can be enabled (might come at the expense of speed)
## see https://software.intel.com/en-us/articles/dapl-ud-support-in-intel-mpi-library
# export I_MPI_DAPL_UD=1

##gather MPI statistics to be analyzed with itac's mps tool
# export I_MPI_STATS=all

##gather MPI debug information (high verbosity)
# export I_MPI_DEBUG=5

# mpiexec -n $NSLOTS octopus_mpi
./start.sh

