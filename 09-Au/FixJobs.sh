#!/bin/bash

set -e

# ListOfJobs="9939236 9939238 9939237 9939233 9939235 9939234 9939220" #Already running
ListOfJobs="9939232 9939219 9939218 9939201 9939217 9939194 9939200 9939195 9939193 9939180 9939190"

echo "$ListOfJobs"

# 1. Search the corresponding paths
for job in ${ListOfJobs}
do
  echo "$job"
  jobpath=$(find -mindepth 2 -maxdepth 3 | grep "$job" | grep 'log' | tail -n1  )
  tail -n2 $jobpath

  jobfolder=$( dirname $jobpath )
  echo $jobfolder

  # Problem was in .bashrc
  cd $jobfolder
  # qsub submit_salomon.pbs
  ${QuantumLaPsources}/09-Au/plotAll.sh
  echo "$job"
  display td.general/Energy.eps
  cd - > /dev/null
done





