#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./octopus.out.%j
#SBATCH -e ./octopus.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J Mission
# Queue (Partition): #general, mpsd (3 days), express, or short
#SBATCH --partition=prace
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=NumNodes
#SBATCH --ntasks-per-node=24
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

pwd >> $QuantumLaPruns/LastRuns.log

#module load intel/16.0 mkl/11.3 impi/5.1.3 gsl/2.1 fftw/3.3.4
module purge
# module load icc/2017.7.259-GCC-6.3.0-2.27
# module load plgrid/tools/intel/17.0.5
module load plgrid/tools/intel/18.0.0
# module load ifort/2017.7.259-GCC-6.3.0-2.27
# module load impi/2017.4.239-iccifort-2017.7.259-GCC-6.3.0-2.27
module load plgrid/tools/impi/2018

# module load imkl/2017.4.239-iimpi-2017c 
module load plgrid/libs/mkl/2018.0.0

# module load libxc/4.2.3-intel-2017a
module load plgrid/libs/libxc/4.2.3

# module load GSL/2.5-intel-2017a
module load plgrid/libs/gsl/2.3

# module load FFTW/3.3.6-intel-2017a 
module load plgrid/libs/fftw/3.3.7

# module load HDF5/1.10.2-intel-2017a
module load plgrid/libs/hdf5/1.10.4

# module load netCDF/4.4.1.1-intel-2017a
module load plgrid/libs/netcdf/4.7.0

module list

# Currently Loaded Modules without libxc
#  1) plgrid/libs/mkl/2018.0.0   5) plgrid/tools/intel/18.0.0
#  2) plgrid/libs/gsl/2.4        6) plgrid/tools/impi/2018
#  3) plgrid/libs/fftw/3.3.7     7) plgrid/libs/hdf5/1.10.4
#  4) plgrid/tools/gcc/6.4.0     8) plgrid/libs/netcdf/4.7.0




# export LIBXC_HOME=/net/people/prace/pr2c0043/libxc
export LIBXC_HOME=/net/software/local/libxc/4.2.3

# export I_MPI_DAPL_TRANSLATION_CACHE=0 #fixes the MPI problem with Draco.

# module load octopus/nightly

# Run the program:
# srun $OCT_HOME/bin/octopus &> octopus-td.$SLURM_JOB_ID.log
mpirun /net/people/prace/pr2c0043/octopus.git/bin/octopus &> octopus-td.$SLURM_JOB_ID.log
# srun /u/tjyderri/octopus.git/bin/octopus &> octopus-td.$SLURM_JOB_ID.log
sleep 5
# sed -i.bak 's/FromScratch/#FromScratch/g' inp #Avoid to kill the data while restarting the job
