#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./octopus.out.%j
#SBATCH -e ./octopus.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J backup
# Queue (Partition): #general, mpsd (3 days), express, or short
#SBATCH --partition=prace  
# plgrid-testing
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#
# Request 120 GB of main Memory per node in Units of MB:
#SBATCH --mem=122000
# Email management
#SBATCH --mail-type=ALL
#SBATCH --mail-user=derrien@fzu.cz
#
# Wall clock limit:
#SBATCH --time=2:00:00
# sbatch -d afternotok:$SLURM_JOB_ID submit.sbatch

pwd >> $QuantumLaPruns/LastRuns.log

#module load intel/16.0 mkl/11.3 impi/5.1.3 gsl/2.1 fftw/3.3.4

# . /etc/profile.d/modules.sh

module purge

source $QuantumLaPsources/octopus_slabs/machines.env
CheckComputerName
LoadEnvLibraries

# export I_MPI_DAPL_TRANSLATION_CACHE=0 #fixes the MPI problem with Draco.

# module load octopus/nightly

# Run the program:
# srun $OCT_HOME/bin/octopus &> octopus-td.$SLURM_JOB_ID.log

# source $QuantumLaPsources/03-Bulk3D/plotAll.sh >& plotting.$SLURM_JOB_ID.log
backup_StH_fixed >& backup-StH.$SLURM_JOB_ID.log
backup_SinSquared >& backup-SinSquared.$SLURM_JOB_ID.log
backup_SinFourth >& backup-SinFourth.$SLURM_JOB_ID.log
echo "Backup is finished."

