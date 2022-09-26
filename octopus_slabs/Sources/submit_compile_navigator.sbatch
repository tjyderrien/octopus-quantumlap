#!/bin/bash -l
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J compile
# Queue (Partition): 
#SBATCH --partition=cpu1
#SBATCH -A frecuencia
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#
# Request 120 GB of main Memory per node in Units of MB:
# S BATCH --mem=95000
# Email management
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kristyna.gazdova@hilase.cz
#
# Wall clock limit:
#SBATCH --time=2:00:00
# sbatch -d afternotok:$SLURM_JOB_ID submit.sbatch

pwd >> $QuantumLaPruns/LastRuns.log


module purge
cat /proc/cpuinfo

module list

# export I_MPI_DAPL_TRANSLATION_CACHE=0 #fixes the MPI problem with Draco.

# Run the program:
# srun $OCT_HOME/bin/octopus &> octopus-td.$SLURM_JOB_ID.log
./compile-octopus-navigator.sh >& octopus-install.$SLURM_JOB_ID.log 
# srun /u/tjyderri/octopus.git/bin/octopus &> octopus-td.$SLURM_JOB_ID.log
sleep 5

