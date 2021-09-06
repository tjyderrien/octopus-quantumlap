#!/bin/bash
set -e 

source $QuantumLaPsources/octopus_slabs/machines.env
CheckComputerName
LoadEnvLibraries

echo "[Check] Variable QuantumLaPsources: $QuantumLaPsources"
echo "[Check] Variable OctopusSlabsSources: $OctopusSlabsSources"

echo "Description: Massive job submitter for QuantumLaP."
echo "WARNING:     1 job ~ 19,200 core-hours per job."
echo " ** HANDLE WITH CARE. ** "
echo ""
echo "Usage: ./flood.sh"
echo "** This file uses run_td.sh."
echo "** DRY-RUN: run_td.sh contains a dry-run for safety."

echo "Note: In this file, you can define all laser pulses of the run you may like to submit."

## ========= COMMON PARTS =============
echo "** Info: defining material and methods"
# Numerics: advise, dont change this! It must be equal to the corresponding ground state. TODO: grasp it automatically in groud state folder. 
#VacuumSize="12." #Optimized for 2D-period MoveIons=No in GS.
# VacuumSize="24." #Hopefully sufficient for TD runs
#VacuumSize="18." #TODO: Optimize for 3D-period MoveIons=Yes.
MeshSizeX="0.23" # 0.20"
# MeshSizeZ=$MeshSizeX
#"0.18" #"0.18" #Octopus-recommanded value: 0.23 Bohr ~ 0.43 Angstroms. After convergence test on current for slab: use 0.18 in Z.

# Number of computer nodes
NumNodes=16

# Sample
SampleLayers=1
SampleXY=1
InitialTemperature="0e0"
CrossPolarization=0 #1

InducedFields="0" #0 or 1: avoid, as understanding is non-trivial there. 
UseLinearResponse="0"
TheoryLevel="dft" #Must be filled: dft #independent_particles /!\ independent-particles != TDFreezeHXC 
Functional="LDA" #LDA defines 'lda_x + lda_c_pz_mod'. TB09 defines "mgga_x_tb09 + lda_c_pw". See http://octopus-code.org/doc/7.0/html/vars.php?page=alpha.
FreezeXC=0 #0 or 1. If 1: TDFreezeHXC = yes is added to the inp. 

echo "** Info: defining the pulse shape"

# delay="1e-15 3e-15 5e-15 7e-15 9e-15 10e-15 20e-15 30e-15" #precising the range
# delay="0e-15" 
delay="0e-15" #sets a delay, unless using the SoftTopHat mode then it's pulse duration.

Mission="HHG-conv" #this replaces the mission in submit.sbatch file. 

# ============ SHIFT CURRENT LASER PARAMETERS =========
#PulseShape="SinSquared"
#lambda1="800e-9"
#lambda2="1850e-9 1950e-9 2050e-9 2150e-9"

# ============ MERDJI LASER PARAMETERS =========

# PROBE PULSE
GaugeVectorFieldAmplitude="0.0" #0: no probe pulse. >0: probe pulse strength in Hartree. 
GaugeFieldDelayValue="0e-15" #delay of probe pulse (in number of cycles?)

PulseShape="SinSquared"
RiseCycle="2" # SYMBOLIC, it is for STH pulses, not SinSquared

lambda1="2656e-9" #not normalized (LDA for 2 um: 2656 nm. TB09 for 2 um: 2237 nm 
lambda2="400e-9" #detuned 2.18: 1422e-9
# lambda2="1030e-9 1050e-9 1070e-9 1090e-9 1110e-9 1130e-9 1150e-9 1170e-9 1190e-9 1210e-9 1230e-9 1250e-9 1270e-9 1290e-9 1310e-9 1330e-9 1350e-9 1370e-9 1390e-9 1410e-9 1430e-9 1450e-9 1470e-9 1490e-9 1510e-9 1530e-9" #MERDJI #seq 1030 20 1550 | awk '{ print $1."e-9" }' | tr "\n" " "

## Field amplitude
Efield1="3E9 9E9" #Exp. 3GV/m (Kozak) 
Efield2="0E0" #"9e9"

## Pulse duration
tau1="13.945e-15"
tau2="30e-15"

## Angle of incidence [for computer cost reasons, we submit 1 ground state per angle, with symmetries. Doing 1 sample for N angles is too costly since symmetries are disabled.]
AngleIncidence_list="0 0.25" 
# AngleIncidence_list="0 0.05 0.10 0.15 0.20 0.25"
# AngleIncidence_list="0.02 0.04 0.06 0.08 0.12 0.14 0.16 0.18 0.22 0.24"
# AngleIncidence_list="0.01 0.02 0.03 0.04 0.06 0.07 0.08 0.09 0.11 0.12 0.13 0.14 0.16 0.17 0.18 0.19 0.21 0.22 0.23 0.24"

# Phase
CEP2="0.0" #"-0.5 -0.25 0.0 0.25 0.5" #0.5 0.333333333 0.25 0.16666667 -0.5 -0.333333333 -0.25 -0.16666667" #will be multiplied by pi in script file
# CEP1="0.5 0.333333333 0.25 0.16666667 0.142857142857143 0.125 0.111111111111111"
CEP1="0.0"

echo "** Info: laser parameters have been defined. Running. "

for iAngleIncidence in ${AngleIncidence_list};
do
for iGaugeVectorFieldAmplitude in ${GaugeVectorFieldAmplitude};
do
for iGaugeFieldDelayValue in ${GaugeFieldDelayValue};
do
for iTheoryLevel in ${TheoryLevel};
do
for iCrossPola in ${CrossPolarization}; 
do
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
#          for iVacuumSize in ${VacuumSize}
#          do
          for iMeshSizeX in ${MeshSizeX}
          do
		iMeshSizeZ=$iMeshSizeX
           for iRiseCycle in ${RiseCycle}
           do
#           for iMeshSizeZ in ${MeshSizeZ}
#           do 
            # Create the working folder
            ## run_td.sh usually does that. 
	    # Prepare the ground state
            # Propagate the pulse in time!
#             echo "Submitting..." 
             $OctopusSlabsSources/03-Bulk3D/run_td.sh ${SampleLayers} ${SampleXY} ${iMeshSizeX} ${iMeshSizeZ} ${iLambda1} ${iLambda2} ${iEfield1} ${iEfield2} ${itau1} ${itau2} ${iCEP1} ${iCEP2} ${iDelay} ${InitialTemperature} ${PulseShape} ${InducedFields} ${iCrossPola} ${iTheoryLevel} ${Functional} ${FreezeXC} ${iRiseCycle} ${Mission} ${NumNodes} ${iGaugeVectorFieldAmplitude} ${iGaugeFieldDelayValue} ${UseLinearResponse} ${iAngleIncidence}
#           done
          done
#          done
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
done
done
done

