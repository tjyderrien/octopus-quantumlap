#!/bin/bash

set -e

source $QuantumLaPsources/octopus_slabs/machines.env
CheckComputerName

echo "Massive job submitter for QuantumLaP."
echo "WARNING: 1 job ~ 19,200 core-hours per job."
echo "HANDLE WITH CARE."
echo "***** MASSIVE GROUND STATE GENERATOR *****"
echo "Usage: ./flood_gs.sh"
echo "This file uses run_gs.sh."
echo "DRY-RUN: run_gs.sh contains a dry-run for safety."

echo "In this file, you can define all numerical features of the run you may like to test."

# Numerics: advise, dont change this! It must be equal to the corresponding ground state. TODO: grasp it automatically in ground state folder.
#VacuumSize="12."
# VacuumSize="24."
# VacuumSize="10."
MeshSizeX_SI="0.20" #"0.15 0.18 0.20 0.23" #0.23 Bohr, octopus-recommended value
#MeshSizeZ=${MeshSizeX_SI} #Convergence study showed that 0.18 is enough to compute the current. #However, it may not be enough to build up the N_exc via using DFT density. 
#MeshSizeZ="0.10"

# Sample
# SampleLayers="6 8 10 12 14 16 18 20"
NumNodes=1
Mission="AuGS"
SampleLayers="1"
SampleXY="1"
kPoints="32" # "12 16 24"
InitialTemperature="0e0"
SmearingTe="0.025" #"1 2 5 10 20" #0.025 eV
TheoryLevel="dft" #"independent_particles" "dft". Empty means dft. Indep_particles != TDFreezeXC !
Functional="PBE" #LDA, PBE, PBEsol, TB09". See http://octopus-code.org/doc/7.0/html/vars.php?page=alpha
Relativistic="no" #yes" #default: no. Optical response for Au fails without relativistic treatment. Spin-Orbit coupling is also important for optical properties. See New Journal of Physics 12 (2010) 103048.

# SampleXY="8 10 12 14 16 18 20"
for iSmearingTe in ${SmearingTe}
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
      $QuantumLaPsources/09-Au/run_gs.sh ${iSampleLayers} ${iSampleXY} ${iMeshSizeX} ${iMeshSizeZ} ${ikPoints} ${InitialTemperature} ${TheoryLevel} ${Functional} ${Mission} ${NumNodes} ${iSmearingTe} ${Relativistic}
 #    done
    done
 #   done
  done
 done
done
done
