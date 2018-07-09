#!/bin/bash

set -e

source $QuantumLaPsources/machines.env
CheckComputerName

echo "Massive job submitter for QuantumLaP."
echo "WARNING: 1 job ~ 19,200 core-hours per job."
echo "HANDLE WITH CARE."
echo "***** MASSIVE GROUND STATE GENERATOR *****"
echo "Usage: ./flood_gs.sh"
echo "This file uses run_gs.sh."
echo "DRY-RUN: run_gs.sh contains a dry-run for safety."

echo "In this file, you can define all numerical features of the run you may like to test."

# Numerics: advise, dont change this! It must be equal to the corresponding ground state. TODO: grasp it automatically in groud state folder. 
#VacuumSize_SI="12."
VacuumSize_SI="24."
# VacuumSize_SI="10."
MeshSizeX_SI="0.23" #Angstroms, octopus-recommended value
MeshSizeZ_SI="0.18" #Convergence study showed that 0.18 is enough to compute the current.
#MeshSizeZ_SI="0.10"

TheoryLevel="dft" #"independent_particles" "dft". Empty means dft. Indep_particles != TDFreezeXC !
Functional="LDA" #LDA or TB09". See http://octopus-code.org/doc/7.0/html/vars.php?page=alpha

# Sample
# SampleLayers="8 10" # 12 14 16 18 20"
SampleLayers="6"
SampleXY="4"
# SampleXY="8 10 12 14 16 18 20"

for iSampleLayers in ${SampleLayers}
do
 for iSampleXY in ${SampleXY}
 do
  for iVacuumSize in ${VacuumSize_SI}
  do
   for iMeshSizeX in ${MeshSizeX_SI}
   do
    for iMeshSizeZ in ${MeshSizeZ_SI}
    do
     # Create the working folder
     ## run_td.sh usually does that. 
     # Prepare the ground state
     # Propagate the pulse in time!
     echo "** SUBMITTING JOB with: SampleLayers=${iSampleLayers} SampleXY=${iSampleXY}" 
     $OctopusSlabsSources/02-Slab2D/run_gs.sh ${iSampleLayers} ${iSampleXY} ${iVacuumSize} ${iMeshSizeX} ${iMeshSizeZ} ${TheoryLevel} ${Functional}
    done
   done
  done
 done
done
