#!/bin/bash
### DESCRIPTION
# Basic example of running a job.
# source spacing.sh >& run.txt

RunningTrunc=/scratch/tjyderri/octopus-svn-15067-nfft/src/main
echo "/!\ Please control version number. octopus_mpi --> octopus for version 15067."
mpiexec -n $NSLOTS $RunningTrunc/octopus >& run.txt
