#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./octopus.out.%j
#SBATCH -e ./octopus.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J QuantumLaP
# Queue (Partition):
#SBATCH --partition=general
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#
# Request 120 GB of main Memory per node in Units of MB:
#SBATCH --mem=122000
# Email management
#SBATCH --mail-type=ALL
#SBATCH --mail-user=derrien@fzu.cz
#
# Wall clock limit:
#SBATCH --time=24:00:00
# sbatch -d afternotok:$SLURM_JOB_ID submit.sbatch

pwd >> $OctopusSlabsRuns/LastRuns.log

#module load intel/16.0 mkl/11.3 impi/5.1.3 gsl/2.1 fftw/3.3.4
module purge
module load intel/17.0
module load impi/2017.3
module load mkl/2017
module load libxc/3.0.0
module load gsl/2.3
module load fftw/3.3.6
module load netcdf-serial/4.4.1.1
module list

export I_MPI_DAPL_TRANSLATION_CACHE=0 #fixes the MPI problem with Draco.

#module load octopus/nightly

# Run the program:
#srun $OCT_HOME/bin/octopus &> octopus-td.$SLURM_JOB_ID.log
srun /u/tjyderri/octopus.git/bin/octopus &> octopus-td.$SLURM_JOB_ID.log
sleep 5
sed -i.bak 's/FromScratch/#FromScracth/g' inp #Avoid to kill the data while restarting the job
