#!/bin/bash

# DRACO@MPCDF remote visualization SLURM submit script.
# 2016, khr@rzg.mpg.de

#SBATCH -o ./job.out.%j
#SBATCH -e ./job.err.%j
#SBATCH -D ./
#SBATCH -J remotevis
#SBATCH --partition=viz
#SBATCH --constraint="gpu"
# --- do not change the following four parameters
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=127000
# ---
#SBATCH --mail-type=none
# --- session duration, up to 24h2#SBATCH --time=01:00:00
#SBATCH --time=10:00:00
module load remotevis

srun ${REMOTEVIS_HOME}/slurm/vis_session

