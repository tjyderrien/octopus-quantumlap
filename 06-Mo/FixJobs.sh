#!/bin/bash

set -e

source ${QuantumLaPsources}/machines.env
CheckComputerName
LoadEnvLibraries

ListOfJobs="9939231 9939146 9939144 9939145 9939142 9939141 9939137 9939138 9939140 9939135 9939134" 

echo "$ListOfJobs"

# 1. Search the corresponding paths
for job in ${ListOfJobs}
do
  echo "$job"
  jobpath=$(find . -mindepth 2 -maxdepth 3 | grep "$job" | grep 'log' | tail -n1  )
  tail -n2 $jobpath

  jobfolder=$( dirname $jobpath )
  echo $jobfolder

  # Problem was in .bashrc
  cd $jobfolder
  echo $jobfolder
  # qsub submit_salomon.pbs
  ${QuantumLaPsources}/09-Au/plotAll.sh
  display td.general/Energy.eps 
  cd - > /dev/null
done





