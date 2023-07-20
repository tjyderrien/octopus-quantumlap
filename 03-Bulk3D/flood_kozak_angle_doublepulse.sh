#!/bin/bash

source $QuantumLaPsources/octopus_slabs/machines.env
CheckComputerName
LoadEnvLibraries

echo "[Check] Variable QuantumLaPsources: $QuantumLaPsources"
echo "[Check] Variable OctopusSlabsSources: $OctopusSlabsSources"

echo "** Info: defining material and methods"
# Numerics: advise, dont change this! It must be equal to the corresponding ground state. TODO: grasp it automatically in groud state folder. 
#VacuumSize="12." #Optimized for 2D-period MoveIons=No in GS.
# VacuumSize="24." #Hopefully sufficient for TD runs
#VacuumSize="18." #TODO: Optimize for 3D-period MoveIons=Yes.
MeshSizeX="0.23" # 0.20"
# MeshSizeZ=$MeshSizeX
#"0.18" #"0.18" #Octopus-recommanded value: 0.23 Bohr ~ 0.43 Angstroms. After convergence test on current for slab: use 0.18 in Z.

# Number of computer nodes

NumNodes=4

# Sample
SampleLayers=1
SampleXY=1
InitialTemperature="0e0"
CrossPolarization=0 #1

InducedFields="0" #0 or 1: avoid, as understanding is non-trivial there. 
UseLinearResponse="0"
TheoryLevel="dft" #Must be filled: dft #independent_particles /!\ independent-particles != TDFreezeHXC 
Functional="TB09" #LDA defines 'lda_x + lda_c_pz_mod'. TB09 defines "mgga_x_tb09 + lda_c_pw". See http://octopus-code.org/doc/7.0/html/vars.php?page=alpha.
FreezeXC=0 #0 or 1. If 1: TDFreezeHXC = yes is added to the inp. 
TimeStep=0.25
echo "** Info: defining the pulse shape"

# delay="1e-15 3e-15 5e-15 7e-15 9e-15 10e-15 20e-15 30e-15" #precising the range
# delay="0e-15" 
delay="0e-15" #sets a delay, unless using the SoftTopHat mode then it's pulse duration.

Mission="HHG-3omega" #this replaces the mission in submit.sbatch file. 

# ============ SHIFT CURRENT LASER PARAMETERS =========
#PulseShape="SinSquared"
#lambda1="800e-9"
#lambda2="1850e-9 1950e-9 2050e-9 2150e-9"

# ============ MERDJI LASER PARAMETERS =========

# PROBE PULSE
GaugeVectorFieldAmplitude="0.0" #0: no probe pulse. >0: probe pulse strength in Hartree. 
GaugeFieldDelayValue="0e-15" #delay of probe pulse (in number of cycles?)

PulseShape="SinFourth"
RiseCycle="2" # SYMBOLIC, it is for STH pulses, not SinSquared

lambda1="2237e-9" #normalized (LDA for 2 um: 2656 nm. TB09 for 2 um: 2237 nm 
# lambda1="407.8e-9" #prepulse for diret TB09 band gap (no need for normalization, as it is already for bandd gap of TB09)
# lambda1="1137.47e-9" #prepulse for indirect TB09 band gap
# lambda1="2656e-9" #normalized (2000 nm, TB09) 
# lambda1="8.388157894736841e-07" #normalized 750 nm, TB09
lambda2="745.66666e-9" #3 omega
# lambda2="1030e-9 1050e-9 1070e-9 1090e-9 1110e-9 1130e-9 1150e-9 1170e-9 1190e-9 1210e-9 1230e-9 1250e-9 1270e-9 1290e-9 1310e-9 1330e-9 1350e-9 1370e-9 1390e-9 1410e-9 1430e-9 1450e-9 1470e-9 1490e-9 1510e-9 1530e-9" #MERDJI #seq 1030 20 1550 | awk '{ print $1."e-9" }' | tr "\n" " "

## Field amplitude


Efield1="1.6E9" # Martin asked for increasing to 1.6E9, from 1.4E9. #"9e9"
# Efield2="0E0"
Efield2="7E7" #"1.4E9" #3E9 9E9" #1.4E9 #I/n this time #3E9" #9E9" #Exp. 3GV/m (Kozak) 
## Pulse duration
# tau1="13.945e-15" #tau1=10 fs with LDA gap renormalized to exp. 3.4 eV gap at 0K. 
# tau1="11.1842e-15" #tau1=10 fs with TB09 gap renormalized to exp. 3.4 eV at 0K.
# tau1="50e-15"
tau1="33.55e-15"
#tau1="33.55e-15" #67.1e-15
tau2=$tau1 #"67.1e-15" #-33.55e-15" #tau=25 fs with TB09 band gap renormalized to 3.4 eV at 0K.
# tau2="11.18e-15" #use if we keep same laser for 1 and 3 omega (normalized to TB09)

# tau1="6.710526315789472e-15" #tau=6 fs, 750 nm, normalized to TB09
# tau1="33.203125e-15" #tau1=25 fs with LDA gap renormalized to exp. 3.4 eV at 0K. 
# tau1="66.40625e-15" #33.203125e-15" #tau1=25 fs with LDA gap renormalized to exp. 3.4 eV at 0K. 
# tau2="30e-15"

## Angle of incidence [for computer cost reasons, we submit 1 ground state per angle, with symmetries. Doing 1 sample for N angles is too costly since symmetries are disabled.]

# SmearingTe="2 5 10 20" #1 (done)
EnableSmearing=0
ExtraStates=0
SmearingTe="0" #1.09 3.06" #"0.025 0.05 0.1 0.5" WE NEED TO FIND A WAY TO DISABLE THIS

# AngleIncidence_list="0.02 0.04 0.06 0.08 0.12 0.14 0.16 0.18 0.22 0.24"
AngleIncidence_list="0" #0.25" # 0" 
# AngleIncidence_list="0" # 0.25" 
# AngleIncidence_list="0.05 0.10 0.15 0.20"
# AngleIncidence_list="0.12"
# AngleIncidence_list="0.04 0.06 0.08 0.12 0.14 0.16 0.18 0.22 0.24"

# AngleIncidence_list="0.01 0.02 0.03 0.04 0.06 0.07 0.08 0.09 0.11 0.12 0.13 0.14 0.16 0.17 0.18 0.19 0.21 0.22 0.23 0.24"

# Phase
CEP1="-0.08569124231348493" #very good prediction for centering pulse1 using tau=33.55e-15 at 2237 nm.
# CEP2="-0.25707372694045255" #ready for centering phase in 30 fs envelope at 3omega (omega=2000 nm).
# CEP1="-0.05738340147065954" #ready for 10 fs pulse
# CEP2="-0.17215092237698829" #ready for 10 fs pulses
# CEP1="0.3064238931585048" #ready for 50 fs pulse
# CEP2="0.9192716794755077" #ready for 50 fs pulses
# CEP1="-0.5 -0.25 0.0 0.25 0.5" #0.5 0.333333333 0.25 0.16666667 -0.5 -0.333333333 -0.25 -0.16666667" #will be multiplied by pi in script file
# CEP1="0.5 0.333333333 0.25 0.16666667 0.142857142857143 0.125 0.111111111111111"
# CEP1="-0.375" #set for 30 fs #No. Error based on A(t), not E(t). 
# CEP1="0.21" #set for 30 fs symmetric for E1(t)
# CEP2="0.0"


CEP2="-1.25707373 -1.15707373 -1.05707373 -0.95707373 -0.85707373 -0.75707373 -0.65707373 -0.55707373 -0.45707373 -0.35707373 -0.25707373 -0.15707373 -0.05707373 0.04292627 0.14292627 0.24292627 0.34292627 0.44292627 0.54292627 0.64292627 0.74292627" #this is -1 pi to 1 pi, but centered exactly.
# CEP1="-1.0 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0"
# CEP1="-0.09 -0.08 -0.07 -0.06 -0.05 -0.04 -0.03 -0.02 -0.01 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09"
# CEP1="-0.69 -0.68 -0.67 -0.66 -0.65 -0.64 -0.63 -0.62 -0.61 -0.59 -0.58 -0.57 -0.56 -0.55 -0.54 -0.53 -0.52 -0.51"
# CEP1="0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29"
if [ ${EnableSmearing} == 0 ]; then
  SmearingTe="0"
fi

echo "** Info: laser parameters have been defined. Running. "
for iSmearingTe in ${SmearingTe};
do
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
             $OctopusSlabsSources/03-Bulk3D/run_td.sh ${SampleLayers} ${SampleXY} ${iMeshSizeX} ${iMeshSizeZ} ${iLambda1} ${iLambda2} ${iEfield1} ${iEfield2} ${itau1} ${itau2} ${iCEP1} ${iCEP2} ${iDelay} ${InitialTemperature} ${PulseShape} ${InducedFields} ${iCrossPola} ${iTheoryLevel} ${Functional} ${FreezeXC} ${iRiseCycle} ${Mission} ${NumNodes} ${iGaugeVectorFieldAmplitude} ${iGaugeFieldDelayValue} ${UseLinearResponse} ${iAngleIncidence} ${EnableSmearing} ${iSmearingTe} ${ExtraStates} ${TimeStep}
#           done
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
done
done

