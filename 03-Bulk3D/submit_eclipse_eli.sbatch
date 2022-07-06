#!/bin/bash  
#SBATCH --account=tderrien-acc1
#SBATCH --partition=queue12h
#SBATCH -N NumNodes # number of nodes 
#SBATCH --clusters=eclipse
#SBATCH -t 0-2:00 # time (D-HH:MM) 
#SBATCH -o octopus.out.%N.%j # STDOUT 
#SBATCH -e octopus.err.%N.%j # STDERR 
#SBATCH --mail-type=END,FAIL # notifications for job done & fail 
#SBATCH --mail-user=derrien@fzu.cz # send-to address

pwd >> $QuantumLaPruns/LastRuns.log

module purge

source ${QuantumLaPsources}/octopus_slabs/modules.env
LoadOctopusEclipseModulesGNU

# module load OpenMPI/1.8.4-GCC-4.8.4
# mpirun -npernode 2  ./stream_mpi_test

octopus_mpiio=""
octopus_latest_gnu="/gpfs/home/tderrien/octopus.gnu.latest/binoctopus" 

mpirun ${octopus_latest_gnu} &> octopus-td.${SLURM_JOB_ID}.log
sleep 5

echo "Done"


