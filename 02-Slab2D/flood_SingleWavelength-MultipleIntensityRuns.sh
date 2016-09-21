#!/bin/bash

echo "Massive job submitter for QuantumLaP."
echo "WARNING: 1 job ~ 19,200 core-hours per job."
echo "HANDLE WITH CARE."

echo "Usage: ./flood.sh"
echo "This file uses run_td.sh."
echo "DRY-RUN: run_td.sh contains a dry-run for safety."

echo "In this file, you can define all laser pulses of the run you may like to submit."
#TODO: develop of version for the ground state, with mesh and vacuum convergence. 

# Numerics: advise, dont change this! It must be equal to the corresponding ground state. TODO: grasp it automatically in groud state folder. 
VacuumSize="10."
MeshSize="0.23"

# Sample
SampleLayers=6
SampleXY=4

# Laser parameters
lambda1="800e-9"
lambda2="400e-9"
Efield1="1E8 5E8 1E9 5E9 1E10 5E10 1E11 5E11 1E12"
Efield2="0E0"
tau1="10e-15"
tau2="10e-15"
CEP1="0.0"
CEP2="0.0"
delay="0.0"

for iLambda1 in ${lambda1};
do
 for iLambda2 in ${lambda2};
 do
  for iEfield1 in ${Efield1}; 
  do
   for iEfield2 in ${Efield2};
   do
    for itau1 in ${tau1};
    do
     for itau2 in ${tau2};
     do
      for iCEP1 in ${CEP1};
      do
       for iCEP2 in ${CEP2};
       do
        for iDelay in ${delay}
        do
         for iVacuumSize in ${VacuumSize}
         do
          for iMeshSize in ${MeshSize}
          do
           # Create the working folder
           ## run_td.sh usually does that. 
	   # Prepare the ground state
           # Propagate the pulse in time! 
           ./run_td.sh ${SampleLayers} ${SampleXY} ${iVacuumSize} ${iMeshSize} ${iLambda1} ${iLambda2} ${iEfield1} ${iEfield2} ${itau1} ${itau2} ${iCEP1} ${iCEP2} ${iDelay}
          done
         done
        done
       done
      done
     done
    done
   done
  done
 done
done
