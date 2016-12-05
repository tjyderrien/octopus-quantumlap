#!/bin/bash

echo "Description: Massive job submitter for QuantumLaP."
echo "WARNING:     1 job ~ 19,200 core-hours per job."
echo " ** HANDLE WITH CARE. ** "
echo ""
echo "Usage: ./flood.sh"
echo "** This file uses run_td.sh."
echo "** DRY-RUN: run_td.sh contains a dry-run for safety."

echo "Note: In this file, you can define all laser pulses of the run you may like to submit."
#TODO: develop of version for the ground state, with mesh and vacuum convergence. 

# Numerics: advise, dont change this! It must be equal to the corresponding ground state. TODO: grasp it automatically in groud state folder. 
#VacuumSize="12." #Optimized for 2D-period MoveIons=No in GS.
VacuumSize="24." #Hopefully sufficient for TD runs
#VacuumSize="18." #TODO: Optimize for 3D-period MoveIons=Yes.
MeshSizeX="0.23"
MeshSizeZ="0.23" #"0.18" #Octopus-recommanded value: 0.23 Bohr ~ 0.43 Angstroms. After convergence test on current for slab: use 0.18 in Z.

# Sample
SampleLayers=1
SampleXY=1

# Laser parameters
lambda1="200e-9 300e-9 400e-9 500e-9 600e-9 700e-9 900e-9 1000e-9 1100e-9 1200e-9"
lambda2="400e-9"
# Field amplitude
# Efield1="1E8 5E8 1E9 5E9 9E9 1E10"
# Efield1="9E9" #9E9 is fluence: 0.1 J/cm2, tau:10 fs. #TODO with scripted intensity
Efield2="0E9"
# Efield1="1E8 5E8 1E9 5E9 1E10 5E10 1E11 5E11 1E12"
Efield1="9E9"
# Pulse duration
# tau1="10e-15 20e-15 30e-15"
tau1="10e-15"
tau2="10e-15"
# Phase
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
          for iMeshSizeX in ${MeshSizeX}
          do
           for iMeshSizeZ in ${MeshSizeZ}
           do 
            # Create the working folder
            ## run_td.sh usually does that. 
	    # Prepare the ground state
            # Propagate the pulse in time! 
            ./run_td.sh ${SampleLayers} ${SampleXY} ${iVacuumSize} ${iMeshSizeX} ${iMeshSizeZ} ${iLambda1} ${iLambda2} ${iEfield1} ${iEfield2} ${itau1} ${itau2} ${iCEP1} ${iCEP2} ${iDelay}
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
done
