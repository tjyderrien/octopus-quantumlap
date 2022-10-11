#!/bin/bash


#set -e

source $QuantumLaPsources/octopus_slabs/machines.env
CheckComputerName
LoadEnvLibraries

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
# VacuumSize="24."
# VacuumSize="10."
MeshSizeX_SI="0.23" #very much used for far. #0.15 0.20" #0.23 Bohr, octopus-recommended value
#MeshSizeZ=${MeshSizeX_SI} #Convergence study showed that 0.18 is enough to compute the current. #However, it may not be enough to build up the N_exc via using DFT density. 
#MeshSizeZ="0.10"

# Sample
# SampleLayers="6 8 10 12 14 16 18 20"

NumNodes=4

Mission="Si_GS_Te"

SampleLayers="1"
SampleXY="1"

# kPoints="4" # 16"
kPoints="40" # 16"
# kPoints="24" # 16"

InitialTemperature="0.025 1 2 5 10" #eV
TheoryLevel="dft" #"independent_particles" "dft". Empty means dft. Indep_particles != TDFreezeXC !
Functional="TB09" #LDA or TB09". See http://octopus-code.org/doc/7.0/html/vars.php?page=alpha
# SampleXY="8 10 12 14 16 18 20"

AngleIncidence_list="0" # 0.25" #"0.05 0.10 0.15 0.20 0.25"

# AngleIncidence_list="0.14 0.16 0.18 0.22 0.24" #"0.05 0.10 0.15 0.20 0.25"  #0.05


# AngleIncidence_list="0.02 0.04 0.06 0.08 0.12 0.14 0.16 0.18 0.22 0.24"
# AngleIncidence_list="0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.11 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 0.20 0.21 0.22 0.23 0.24"

for iInitialTemperature in ${InitialTemperature}
do
for iAngleIncidence in ${AngleIncidence_list}
do
for ikPoints in ${kPoints}
do
 for iSampleLayers in ${SampleLayers}
 do
  for iSampleXY in ${SampleXY}
  do
 #   for iVacuumSize in ${VacuumSize}
 #   do
    for iMeshSizeX in ${MeshSizeX_SI}
    do
     iMeshSizeZ=$iMeshSizeX
 #    for iMeshSizeZ in ${MeshSizeZ}
 #    do
      # Create the working folder
      ## run_td.sh usually does that. 
      # Prepare the ground state
      # Propagate the pulse in time!
      echo "** SUBMITTING JOB with: SampleLayers=${iSampleLayers} SampleXY=${iSampleXY}" 
      $QuantumLaPsources/03-Bulk3D/run_gs.sh ${iSampleLayers} ${iSampleXY} ${iMeshSizeX} ${iMeshSizeZ} ${ikPoints} ${iInitialTemperature} ${TheoryLevel} ${Functional} ${Mission} ${NumNodes} ${iAngleIncidence}
 #    done
    done
   done
  done
 done
done
done
