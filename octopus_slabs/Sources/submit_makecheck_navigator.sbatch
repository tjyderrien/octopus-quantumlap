#!/bin/bash

#SBATCH -A frecuencia
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J make-check
# Queue (Partition): #general, mpsd (3 days), express, or short
#SBATCH --partition=cpu1
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#
# Request 120 GB of main Memory per node in Units of MB:
# S BATCH --mem=122000
# Email management
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kristyna.gazdova@hilase.cz
#
# Wall clock limit:
#SBATCH --time=2:00:00
# s batch -d afternotok:$SLURM_JOB_ID submit.sbatch
# Standard output and error:
# S BATCH -o ./octopus.out.%j
# S BATCH -e ./octopus.err.%j

module purge

module list

module load iccifort/2019.5.281

module load iimpi/2019b

source /apps/software/imkl/2019.5.281-iimpi-2019b/bin/compilervars.sh intel64

module load imkl/2019.5.281-iimpi-2019b
export MKL_HOME="/apps/software/imkl/2019.5.281-iimpi-2019b/"

module load libxc/4.3.4-iccifort-2019.5.281
export LIBXC_HOME="/apps/software/libxc/4.3.4-iccifort-2019.5.281"

module load GSL/2.6-iccifort-2019.5.281
export GSL_HOME="/apps/software/GSL/2.6-iccifort-2019.5.281/"


# export LIBXC_HOME=/net/people/prace/pr2c0043/libxc
# export LIBXC_HOME=/net/software/local/libxc/4.2.3

module load libyaml/0.2.2-GCCcore-8.3.0
export YAML_HOME=/apps/software/libyaml/0.2.2-GCCcore-8.3.0/

# export I_MPI_DAPL_TRANSLATION_CACHE=0 #fixes the MPI problem with Draco.

# Run the program:
octopus_latest="/veracruz/home/p/pr2c0043/octopus.intel/bin/octopus"
octopus_mpiio="/veracruz/home/p/pr2c0043/octopus.mpiio.July/bin/octopus"
octopus_barbora="/veracruz/home/p/pr2c0043/octopus.barbora/bin/octopus"
# srun $OCT_HOME/bin/octopus &> octopus-td.$SLURM_JOB_ID.log
#which mpirun

ldd $octopus_barbora

export OCT_TEST_NJOBS=8
make check #mpirun ${octopus_barbora} &> octopus-td.${SLURM_JOB_ID}.log

# srun /u/tjyderri/octopus.git/bin/octopus &> octopus-td.$SLURM_JOB_ID.log
sleep 5
# sed -i.bak 's/FromScratch/#FromScratch/g' inp #Avoid to kill the data while restarting the job
