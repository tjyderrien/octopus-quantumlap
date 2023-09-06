#!/bin/bash  
#SBATCH --account=hilase
#SBATCH --partition=T1H
#SBATCH -N NumNodes # number of nodes
#SBATCH -n 24 #number of threads per node
#SBATCH --job-name=Mission
# S BATCH --clusters=sunrise
#SBATCH -t 0-01:00 # time (D-HH:MM) 
#SBATCH -o octopus.out.%N.%j # STDOUT 
#SBATCH -e octopus.err.%N.%j # STDERR 
#SBATCH --mail-type=END,FAIL # notifications for job done & fail 
#SBATCH --mail-user=derrien@fzu.cz # send-to address

pwd >> $QuantumLaPruns/LastRuns.log

module purge

source ${QuantumLaPsources}/octopus_slabs/modules.env
LoadOctopusModulesSunrise
# Let's try it manually

# octopus_mpiio="/gpfs/home/tderrien/octopus.parallel_restart.gcc7.libxc4/octopus/_build/src/octopus"
octopus_mpiio="/gpfs/home/tderrien/octopus.parallel_restart.gnu7.libxc4.mpiio/bin/octopus"
# octopus_latest_gnu="/gpfs/home/tderrien/octopus.gnu.latest/bin/octopus" 

srun ${octopus_mpiio} &> octopus-td.${SLURM_JOB_ID}.log
sleep 5

echo "Done"


