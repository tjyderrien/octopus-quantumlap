#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./octopus.out.%j
#SBATCH -e ./octopus.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J QuantumLaP-DBbuild
# Queue (Partition): #general, mpsd (3 days), express, or short
#SBATCH --partition=plgrid-testing
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#
# Request 120 GB of main Memory per node in Units of MB:
#SBATCH --mem=122000
# Email management
#SBATCH --mail-type=ALL
#SBATCH --mail-user=derrien@fzu.cz
#
# Wall clock limit:
#SBATCH --time=1:00:00
# sbatch -d afternotok:$SLURM_JOB_ID submit.sbatch

sleep 5

pwd >> $QuantumLaPruns/LastRuns.log

#module load intel/16.0 mkl/11.3 impi/5.1.3 gsl/2.1 fftw/3.3.4

# . /etc/profile.d/modules.sh

# TEMPORARY FIX
module purge

source $QuantumLaPsources/octopus_slabs/machines.env
CheckComputerName
# export I_MPI_DAPL_TRANSLATION_CACHE=0 #fixes the MPI problem with Draco.

# module load octopus/nightly

# Run the program:
# srun $OCT_HOME/bin/octopus &> octopus-td.$SLURM_JOB_ID.log

# source $QuantumLaPsources/03-Bulk3D/plotAll.sh >& plotting.$SLURM_JOB_ID.log

source ${QuantumLaPsources}/Cleaning/CountInodes.sbatch >& inodes.log 

# ${QuantumLaPsources}/02-Slab2D/plotFinalQuantities.py >& databasing.${SLURM_JOB_ID}.log

