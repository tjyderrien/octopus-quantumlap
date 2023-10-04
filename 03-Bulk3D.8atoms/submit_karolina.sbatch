#!/bin/bash
#SBATCH --partition=qcpu
#SBATCH -J Mission
#SBATCH --mail-user=USER_EMAIL
#SBATCH --mail-type=ALL
#SBATCH -o ./octopus.out.%j
#SBATCH -e ./octopus.err.%j
#SBATCH -A OPEN-29-56
# P B S -l select=NumNodes:ncpus=128:mpiprocs=128:ompthreads=1
#SBATCH --nodes=NumNodes
#SBATCH --ntasks-per-node=1
#SBATCH --time=2-00:00:00

# cd $PBS_O_WORKDIR

pwd >> $QuantumLaPruns/LastRuns.log

source ${QuantumLaPsources}/octopus_slabs/modules.env
LoadOctopusKarolinaIntel
# LD_LIBRARY_PATH=/scratch/tjyderriien/nfft/3.2.4/lib:$LD_LIBRARY_PATH

# module load octopus/5.0.0

## to save memory at high core counts, the "connectionless user datagram protocol" 
## can be enabled (might come at the expense of speed)
## see https://software.intel.com/en-us/articles/dapl-ud-support-in-intel-mpi-library
# export I_MPI_DAPL_UD=1

# export I_MPI_DAPL_TRANSLATION_CACHE=0 #fixes an MPI problem with Draco.

##gather MPI statistics to be analyzed with itac's mps tool
# export I_MPI_STATS=all

##gather MPI debug information (high verbosity)
# export I_MPI_DEBUG=5

# mpiexec -n $NSLOTS octopus_mpi
# ./start.sh
# RunningTrunc=/scratch/tjyderri/octopus-svn-15067/src/main

# JUST TO EXECUTE ONE RUN
MPIIOtrunc=/home/${LOGNAME}/.local/easybuild/software/Octopus/parallel_restart-intel-2020b-mpi/bin/
# srun $OCT_HOME/bin/octopus &> octopus-td.$SLURM_JOB_ID.log
# srun /u/tjyderri/octopus.git/bin/octopus &> octopus-td.$SLURM_JOB_ID.log
# sleep 5
# sed -i.bak 's/FromScratch/#FromScratch/g' inp #Avoid to kill the data while restarting the job
# ~                                                           
mpirun $MPIIOtrunc/octopus >& octopus-td.${SLURM_JOB_ID}.log #If $NSLOTS fails, use PBS_NUM_PPN
