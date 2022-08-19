#!/bin/bash
# Standard output and error:
#SBATCH -o ./octopus.out.%j
#SBATCH -e ./octopus.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J Mission
# Queue (Partition): #general, mpsd (3 days), express, or short
#SBATCH --partition=cpu1
#SBATCH -A frecuencia
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=NumNodes
#SBATCH --ntasks-per-node=48
#
# Request 120 GB of main Memory per node in Units of MB:
# S BATCH --mem=122000
# Email management
#SBATCH --mail-type=ALL
#SBATCH --mail-user=derrien@fzu.cz
#
# Wall clock limit:
#SBATCH --time=48:00:00
# sbatch -d afternotok:$SLURM_JOB_ID submit.sbatch

pwd >> $QuantumLaPruns/LastRuns.log

module purge

source ${QuantumLaPsources}/octopus_slabs/modules.env
LoadOctopusNavigatorModulesFortran2021

# export I_MPI_DAPL_TRANSLATION_CACHE=0 #fixes the MPI problem with Draco.
# https://community.intel.com/t5/Intel-oneAPI-HPC-Toolkit/Memory-registration-cache-feature-in-DAPL-gt-random-failure-in/td-p/1135455

# Run the program:
octopus_latest="/veracruz/home/p/pr2c0043/octopus.intel/bin/octopus"
octopus_mpiio="/veracruz/home/p/pr2c0043/octopus.mpiio.July/bin/octopus"
octopus_barbora="/veracruz/home/p/pr2c0043/octopus.barbora/bin/octopus"
# srun $OCT_HOME/bin/octopus &> octopus-td.$SLURM_JOB_ID.log #AVOID THIS. mpirun works very well. 
mpirun ${octopus_mpiio} &> octopus-td.${SLURM_JOB_ID}.log #priorities: srun, then mpirun

# srun /u/tjyderri/octopus.git/bin/octopus &> octopus-td.$SLURM_JOB_ID.log
sleep 5
#bak 's/FromScratch/#FromScratch/g' inp #Avoid to kill the data while restarting the job
