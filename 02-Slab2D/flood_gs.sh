#!/bin/bash

echo "Massive job submitter for QuantumLaP."
echo "WARNING: 1 job ~ 19,200 core-hours per job."
echo "HANDLE WITH CARE."
echo "***** MASSIVE GROUND STATE GENERATOR *****"
echo "Usage: ./flood_gs.sh"
echo "This file uses run_gs.sh."
echo "DRY-RUN: run_gs.sh contains a dry-run for safety."

echo "In this file, you can define all numerical features of the run you may like to test."

# Numerics: advise, dont change this! It must be equal to the corresponding ground state. TODO: grasp it automatically in groud state folder. 
#VacuumSize="12."
VacuumSize="24."
# VacuumSize="10."
MeshSizeX="0.23" #Bohr, octopus-recommended value
MeshSizeZ="0.18" #Convergence study showed that 0.18 is enough to compute the current.
#MeshSizeZ="0.10"

# Sample
# SampleLayers="8 10" # 12 14 16 18 20"
SampleLayers="6"
SampleXY="6"
# SampleXY="8 10 12 14 16 18 20"

for iSampleLayers in ${SampleLayers}
do
 for iSampleXY in ${SampleXY}
 do
  for iVacuumSize in ${VacuumSize}
  do
   for iMeshSizeX in ${MeshSizeX}
   do
    for iMeshSizeZ in ${MeshSizeZ}
    do
     # Create the working folder
     ## run_td.sh usually does that. 
     # Prepare the ground state
     # Propagate the pulse in time!
     echo "** SUBMITTING JOB with: SampleLayers=${iSampleLayers} SampleXY=${iSampleXY}" 
     ./run_gs.sh ${iSampleLayers} ${iSampleXY} ${iVacuumSize} ${iMeshSizeX} ${iMeshSizeZ}
    done
   done
  done
 done
done
