#!/bin/bash
set -e 

source $QuantumLaPsources/machines.env
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
# Numerics: advise, dont change this! It must be equal to the corresponding ground state. TODO: grasp it automatically in groud state folder. 
#VacuumSize="12." #Optimized for 2D-period MoveIons=No in GS.
# VacuumSize="24." #Hopefully sufficient for TD runs
#VacuumSize="18." #TODO: Optimize for 3D-period MoveIons=Yes.
MeshSizeX="0.23" #0.23" #In general: 0.23 was used."
# MeshSizeZ=$MeshSizeX
#"0.18" #"0.18" #Octopus-recommanded value: 0.23 Bohr ~ 0.43 Angstroms. After convergence test on current for slab: use 0.18 in Z.

# Number of NumNodes
NumNodes=8

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

### MISSION #####
# Mission="Drude"

# Mission="QLaP-1600"
# Mission="3200nm-LinearResponse"
# Mission="CEP-tunneling"
# Mission="InducedField"
# Mission="MidIRGapModification"
# Mission="Si-7PI-3200nm"
# Mission="Yamada0.8eV"
Mission="FigS2-20fs"
# Mission="Fig2-20fsSat"


# Phase
CEP1="0.0" #-0.5 -0.25 0.0 0.25 0.5 0.333333333 0.16666667 -0.333333333 -0.16666667" #will be multiplied by pi in script file
CEP2="0.0" #0.0 0.5 0.333333333 0.25 0.16666667 0.142857142857143 0.125"
AngleOfIncidence="0"
# delay="1e-15 3e-15 5e-15 7e-15 9e-15 10e-15 20e-15 30e-15" #precising the range
# delay="0e-15" 
delay="0.25e-15" #sets a delay, unless using the SoftTopHat mode then it's pulse duration.

# PROBE PULSE
GaugeVectorFieldAmplitude="0.0" #0: no probe pulse. >0: probe pulse strength in Hartree. 
GaugeFieldDelayValue="0e-15" #delay of probe pulse (in number of cycles?)

## ========== ZHUKOV LASER PARAMETERS ============
PulseShape="StH" 
# StH: corrected one. 
# StH4: soft top hat with 4th order sinus on ramp up/down
# SinFourth: sin4(t)
# STH: old buggy STH. 
# SinSquared: sin2(t) function. 
#Zhukov: STH, Merdji: SinSquared. #Gaussian #SinSquared #TopHat #STH = SoftTopHat: WARNING this pulse is defined by the method tdf_from_file. Other methods are generated in octopus directly. 
# SinPower=4
RiseCycle="2" # "1 2 3 4 5 9" #elongates the pulse duration at some point
# lambda1="800e-9"  # "900e-9" # 1100e-9 1300e-9 1400e-9 1500e-9 1600e-9 1700e-9 1800e-9 1900e-9 2000e-9 2100e-9" #BICOLOR SEQUENCE
# lambda1="100e-9 150e-9 200e-9 250e-9 300e-9 350e-9 400e-9 450e-9 500e-9 550e-9 600e-9 650e-9 700e-9 750e-9 850e-9 900e-9 1000e-9 1100e-9 1200e-9 1300e-9 1400e-9 1500e-9 1600e-9 1700e-9 1800e-9 1900e-9 2000e-9"
# lambda2="800e-9"
# lambda2="515e-9 700e-9 900e-9 1030e-9 1060e-9 1100e-9 1300e-9 1500e-9 1700e-9 1900e-9 2100e-9"
# lambda2="515e-9 600e-9 700e-9 800e-9 900e-9 1000e-9 1030e-9 1060e-9 1100e-9 1200e-9 1300e-9 1400e-9 1500e-9 1600e-9 1700e-9 1800e-9 1900e-9 2000e-9 2100e-9" #ZHUKOV

# Field amplitude
Efield_fast="5E8" # 1E9 5E9 1E10"
# Efield_detailed="9.75E8 1E9 1.5E9 2E9"
Efield_detailed="2E8 2.5E8 3E8 3.5E8 4E8 4.5E8" # 5.5E8 6E8 6.5E8 7E8 7.5E8 8E8 8.5E8 9E8 9.5E8 2.5E9 3.5E9 4.5E9 5.5E9 6E9 6.5E9 7E9 7.5E9 8E9 8.5E9 9E9 9.5E9" # is fluence: 0.1 J/cm2, tau:10 fs. #TODO with scripted intensity
Efield_superdetailed="2.1E9 2.2E9 2.3E9 2.4E9 2.5E9 2.6E9 2.7E9 2.8E9 2.9E9 3.1E9 3.2E9 3.3E9 3.4E9" 
Efield_saturation="1E10 2E10 3E10 4E10 5E10 6E10 7E10 8E10 9E10"
Efield1=${Efield_detailed}
# Efield1="9.95E8 1E9 1.05E9"
# Efield1="5E9 5.5E9 6E9 6.5E9 7E9 7.5E9 8E9 8.5E9 9E9"
# Efield1="$Efield_fast $Efield_detailed"
# Efield1="9.2E9 9.3E9 9.4E9 9.5E9 9.6E9 9.7E9 9.8E9 9.9E9" #1.1E9 1.2E9 1.3E9 1.4E9 1.6E9"
# Efield1="1.5E8 2.5E8 3.5E8 4.5E8 5.5E8 6.5E8 7.5E8 8.5E8 9.5E8 1.5E9 2.5E9 3.5E9 4.5E9 5.5E9 6.5E9 7.5E9 8.5E9 9.5E9"
# Efield1="1E10 3E10 4E10 5E10 6E10 7E10 8E10 9E10"
# Efield1="1E11 2E11 3E11 4E11 5E11 6E11 7E11 8E11 9E11"
# Efield1="5.5E8 5E8 6E8 6.5E8 7E8 7.5E8 8.5E8 8E8 9E8 9.5E8 1.5E9 2E9 2.5E9 3E9 3.5E9 4E9 4.5E9 6E9 7E9 8E9 9E9" #Exceptionally, run a series of runs which did not finish. 
# Efield1="9.5E9 1.5E10 2.5E10 3.5E10 4.5E10 5.5E10 6.5E10 7.5E10 8.5E10 9.5E10"
# Efield1="3E8 4E8 5.5E9 6.5E9 7.5E9 8.5E9 9.5E9" 
# Efield1="1E7 5E7 1E8 2E8 3E8" # Super low intensity runs, to be better converged.

# Series of runs for PRL Fig3
# Efield1="5E8 7.5E9 1E9 2E9 3E9"
# Efield1="5E8 3E9 4E9"
# Efield1="5E10 1E11 5E11"
# Efield1="1E8 2.5E8" # 3E9 4E9"
# Efield1="0E0" #"3E9" #2.5E8" # 3E9 4E9"

# STARK SHIFT VALUES.
# Efield1="5E8 7.5E8 1E9 1.5E9 2E9 3E9 4E9 5E9" 
# Efield2="0E0 1E8 5E8 1E9 5E9 1E10 5E10"
# Efield2="0E0"
# tau1="20E-15"
# tau2="20E-15"

## PRL (intensity) figure
# Checking mid-IR band gap at 5-10 um
# tau1="60E-15" # 20e-15 30E-15 10E-15"
# tau1="10.674E-15" #(1 optical cycle at 3200 nm)
tau1="20E-15"
# tau1="10e-15 20e-15 30e-15 40e-15 50e-15"
tau2="20E-15"
# tau2="10.674E-15"
Efield2="0E0"

# lambda1="1600e-9"
# lambda1="200e-9 300e-9 343e-9 400e-9 500e-9 515e-9 600e-9 700e-9 800e-9 900e-9 1000e-9 1030e-9 1064e-9 1100e-9 1200e-9 1300e-9 1400e-9 1500e-9 1600e-9 1700e-9 1800e-9 1900e-9 2000e-9 2100e-9 2600e-9" # 3200e-9"
# lambda1="200e-9 300e-9 343e-9 400e-9 500e-9 515e-9 600e-9 700e-9 800e-9 900e-9 1000e-9 1030e-9 1064e-9 1100e-9 1200e-9 1300e-9 1400e-9 1500e-9 1600e-9 1700e-9 1800e-9 1900e-9 2000e-9 2100e-9"
# energies="0.50 0.54 1.27 1.28 1.29 1.30 1.35 1.90 2.17 2.77 3.83"
# lambda1="2470e-9 3100e-9 4100e-9 6200e-9 12400e-9"
# lambda1="323e-9 448e-9 571e-9 652e-9 918e-9 954e-9 961e-9 968e-9 976e-9 2296e-9 2480e-9"
# lambda1="550e-9 750e-9 850e-9 1150e-9 1250e-9 1350e-9 1450e-9 1650e-9 1750e-9" 
# lambda1="483e-9 484e-9" #these are harmonics of the LDA band gap. 

# lambda1="2200e-9 2600e-9 3200e-9 3600e-9"
# lambda1="1600e-9 2200e-9 3200e-9"
# lambda1="1540e-9" #3200e-9" #Yabana wavelength for linear interband: 1.54 um. 
# lambda1="6.4e-6 8e-6 10e-6 12.4e-6 16e-6 19.2e-6 22.4e-6 25e-6" #how to see the band gap dynamics? See the photons between 5 and 10 um.
# lambda1="1600e-9 2200e-9 3200e-9"
# lambda1="990e-9 1010e-9 1020e-9 1030e-9 1040e-9 1050e-9 1060e-9 1070e-9 1080e-9 1090e-9"
# lambda1="407.58e-9" #1st harmonic of the band gap LDA: 483.70 nm. TB09: 407.58 nm.
lambda1="1600e-9"
lambda2="400e-9"
# lambda2="1030e-9 1600e-9"

# ============ SHIFT CURRENT LASER PARAMETERS =========
#PulseShape="SinSquared"
#lambda1="800e-9"
#lambda2="1850e-9 1950e-9 2050e-9 2150e-9"

# ============ MERDJI LASER PARAMETERS =========
# PulseShape="SinSquared" 

# lambda1="3100e-9" #HAMED wavelength
# lambda2="1550e-9" #detuned 2.18: 1422e-9
# lambda2="1030e-9 1050e-9 1070e-9 1090e-9 1110e-9 1130e-9 1150e-9 1170e-9 1190e-9 1210e-9 1230e-9 1250e-9 1270e-9 1290e-9 1310e-9 1330e-9 1350e-9 1370e-9 1390e-9 1410e-9 1430e-9 1450e-9 1470e-9 1490e-9 1510e-9 1530e-9" #MERDJI #seq 1030 20 1550 | awk '{ print $1."e-9" }' | tr "\n" " "

## Field amplitude
# Efield1="9e9" #MERDJI
# Efield2="9e9"

## Pulse duration
# tau1="30e-15"
# tau2="30e-15"


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
            $OctopusSlabsSources/03-Bulk3D/run_td.sh ${SampleLayers} ${SampleXY} ${iMeshSizeX} ${iMeshSizeZ} ${iLambda1} ${iLambda2} ${iEfield1} ${iEfield2} ${itau1} ${itau2} ${iCEP1} ${iCEP2} ${iDelay} ${InitialTemperature} ${PulseShape} ${InducedFields} ${iCrossPola} ${iTheoryLevel} ${Functional} ${FreezeXC} ${iRiseCycle} ${Mission} ${NumNodes} ${iGaugeVectorFieldAmplitude} ${iGaugeFieldDelayValue} ${UseLinearResponse} ${AngleOfIncidence}
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
