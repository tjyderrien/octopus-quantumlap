#!/bin/bash
# Input: $1: SampleLayers, $2: SampleXY
# Purpose: TD-evolution. 

# set -e
source machines.env
source functions.sh
# ENABLE DRY RUN ? 
DryRun=1
MoveIons=0
kPoints=1

# Renaming inputs to replace into the INP file.
SampleLayers=$1
SampleXY=$2
VacuumSizeInput=$3
MeshSizeXInput=$4
MeshSizeZInput=$5
lambda1=$6
lambda2=$7
Efield1=$8
Efield2=$9
tau1=${10}
tau2=${11}
CEP1=${12}
CEP2=${13}
delay=${14}

ProcessArguments () { # Conversion factors for sample
  Bohr2Angstrom=0.52917721067 #AngstromPerBohr

  # Renaming inputs
  SampleLayers=$1
  SampleXY=$2
  SizeOfVacuumAng=${VacuumSizeInput} #Angstrom
  SizeOfVacuum=$( python_calc "${SizeOfVacuumAng} / ${Bohr2Angstrom}" ) #Bohr
  SizeOfVacuumName=$( python_calc "int(round( ${SizeOfVacuumAng} ))" )
  #MeshSizeAng=0.23 #Angstrom
  MeshSizeXAng=${MeshSizeXInput} #0.23 #Angstrom
  MeshSizeZAng=${MeshSizeZInput} #should be smaller, as we have a surface
  MeshSizeX=$( python_calc "${MeshSizeXAng} / ${Bohr2Angstrom}"  ) #Bohr
  MeshSizeZ=$( python_calc "${MeshSizeZAng} / ${Bohr2Angstrom}"  ) #Bohr
#  MeshSize=$( python_calc "${MeshSizeAng} / ${Bohr2Angstrom}"  ) #Bohr
}

DefineParameters () {
  # Input laser parameters
  RelaxationTime=$( python_calc "0e0"  )

# echo "Generating pulses..."

  echo "TODO: for now, the Python laser generator only makes nice pictures."
#  ../Pulses/BuildPulseDatabase.py > Pulse.dat
#  lambda1=$( cat 'Pulse.dat' | grep 'Wavelength' | awk '{ print $2 }' ) 
#  lambda2=$( cat 'Pulse.dat' | grep 'Wavelength' | awk '{ print $5 }' ) 
#  Efield1=$( cat 'Pulse.dat' | grep 'Amplitude' | awk '{ print $2 }' ) 
#  Efield2=$( cat 'Pulse.dat' | grep 'Amplitude' | awk '{ print $5 }' ) 
#  tau1=$( cat 'Pulse.dat' | grep 'Duration' | awk '{ print $2 }' ) 
#  tau2=$( cat 'Pulse.dat' | grep 'Duration' | awk '{ print $5 }' ) 
#  CEP1=$( cat 'Pulse.dat' | grep 'Phase' | awk '{ print $2*2.*3.141592653/180. }' ) 
#  CEP2=$( cat 'Pulse.dat' | grep 'Phase' | awk '{ print $5*2.*3.141592653/180. }' ) 
#  delay=$( cat 'Pulse.dat' | grep 'Delay' | awk '{ print $5 }' )
  
  sigmaTau=$( awk 'END { printf("%1.15lg\n", 2e0*sqrt(2e0*log(2e0)) )}' /dev/null )
#  TimeAUtoSI=$( python_calc " 1.054571800e-34 / 1.6021766208e-19 / 27.21138505e0 " )
#  TimeStepAU=$( grep 'TDTimeStep' inp-td | awk '{ print $3 }'  )
#  TimeStepSI=$( python_calc " ${TimeStepAU} * ${TimeAUtoSI} " ) 
  ComputingTime=$( python_calc " float( (8.*${tau1}/${sigmaTau} + ${delay} + ${RelaxationTime}) ) " )
  echo "Pulses generated."

  MoveIonsChain=""
  if [ ${MoveIons} == 1 ]; then
	MoveIonsChain="-MoveIons"
  else
        MoveIonsChain=""
  fi
  export groundStatePath=${basePath}/Si-${SampleXY}x${SampleXY}x${SampleLayers}-k-${kPoints}x${kPoints}x1-L${SizeOfVacuumName}-Mesh${MeshSizeX}-Mesh${MeshSizeZ}
  export localPath=${groundStatePath}${MoveIonsChain}/"lambda"${lambda1}-"tau"${tau1}-"field"${Efield1}-"CEP"${CEP1}/"lambda"${lambda2}-"tau"${tau2}-"field"${Efield2}-"CEP"${CEP2}-"delay"${delay}
  if [ ! -d ${localPath} ]; then
	mkdir -p ${localPath}
  fi
}

print_usage () {
  if [ -e ${SampleXY} ]; then
    echo "Usage: "
    echo "./run.sh <Number of Z blocks> <Number of XY atoms> <Vacuum Size (Angstrom), Mesh size XY (Angstrom), Mesh size Z (Ang.)> <lambda1, lambda2 (m)>, <Efield1, Efield2 (V/m)>, <tau1, tau2 (s)>, <CEP1 (rad), CEP2 (rad)>, <delay (s)>"
    exit
  else
    echo "** Sample size: Z=${SampleLayers} layers. XY=${SampleXY} atoms.**"
    echo "** Mesh size: Z=${MeshSizeZInput}. X=${MeshSizeXInput}. **"
  fi
}

ReplaceInputFile () {
  # Define the material parameters
  # We would like to start the td-calculation
  # INP file is modified in root, then copied to the ${localPath} folder  
  if [ -e inp ]; then
	  echo "Input file INP was renamed to .BAK..."
	  mv inp inp.bak
  fi
  echo "Input file INP was set to TIME-DEPENDENT RUN"
  cp inp-td inp
  echo "REPLACING MATERIAL PARAMETERS..."
  sed -i.bak s/SizeOfTheBox/${SampleXY}/g inp
  sed -i.bak s/NumberOfBlocks/${SampleLayers}/g inp
  sed -i.bak s/SizeOfVacuum/${SizeOfVacuum}/g inp
  sed -i.bak s/MeshSizeX/${MeshSizeX}/g inp
  sed -i.bak s/MeshSizeZ/${MeshSizeZ}/g inp
  sed -i.bak s/kPointsValue/${kPoints}/g inp
  echo "REPLACING LASER PULSE PARAMETERS..."
  sed -i.bak s/lambda1Value/${lambda1}/g inp
  sed -i.bak s/lambda2Value/${lambda2}/g inp
  sed -i.bak s/Efield1amplitude/${Efield1}/g inp
  sed -i.bak s/Efield2amplitude/${Efield2}/g inp
  sed -i.bak s/tau1ToReplace/${tau1}/g inp
  sed -i.bak s/tau2ToReplace/${tau2}/g inp
  sed -i.bak s/pulseDelayToReplace/${delay}/g inp
  sed -i.bak s/CEP1replaceMe/${CEP1}/g inp
  sed -i.bak s/CEP2replaceMe/${CEP2}/g inp
  sed -i.bak s/ComputingTimeValueToReplace/${ComputingTime}/g inp

  if [ ${MoveIons} == 1 ]; then
	echo "MoveIons? Addind the function in Input file..."
	echo "MoveIons = yes" >> inp
  fi

}

BuildSample() {
  # Build the sample
  if [ -e Si.xyz ]; then
    rm Si.xyz
  fi
  $AbsPath/Samples/Si100-H/write-geom.py ${SampleLayers} ${SampleXY} | $AbsPath/Samples/fhi2octopus.sh > Si.xyz

  # Add the number of lines at the beginning of the file
  numLines=$( wc -l Si.xyz | awk '{ print $1 }' )

  echo $numLines >> Si.tmp
  echo "# This line is a comment." >> Si.tmp
  cat Si.xyz >> Si.tmp
  cat Si.tmp
  mv Si.tmp Si.xyz
}

SubmitTheJob () {
  echo "== CHECKING GROUND STATE INTEGRITY =="
  #If ground state exists, we shall copy the ground state in the target folder!
  #else: don't run the ground, just stop the program and tell it. 
  if [ -e ${groundStatePath}-Save ]; then
          ./DuplicateGroundState.sh ${groundStatePath}-Save/ ${groundStatePath}${MoveIonsChain}/
  else
    echo "GROUND STATE DOES NOT EXIST. Run it manually using run_gs.sh ${SampleLayers} ${SampleXY}."
    exit
  fi 

  echo "== Copying the file in the submission folder... =="
  echo "Source: ${basePath}."
  echo "Target: ${localPath}."
  cp ${basePath}/inp ${localPath}/
#  cp ${basePath}/submit_mpi.sge ${localPath}/
  cp ${basePath}/Si.xyz ${localPath}/
  cp ${basePath}/$SubmitScript ${localPath}/

  # Start the run
  if [ ${DryRun} == 0 ]; then
          cd ${localPath}
          ${SubmitCommand} ${SubmitScript} > submit.log
	  jobnumber=$( cat submit.log | awk '{ print $4 }') 
	  cp submit.log $jobnumber.log
          cd -
  else
          echo "** DRY-RUN: Change value into the script to unlock the safety. **"
          echo "${ComputerName}: ${SubmitCommand} ${SubmitScript}"
  fi
}

CheckComputerName
ProcessArguments ${SampleLayers} ${SampleXY}
print_usage

# For each calculation to perform
DefineParameters
ReplaceInputFile #this is performed in the common directory!
BuildSample #also performed in directory
# ReplaceInputFile 
SubmitTheJob #then the file are move in the directory where calculation will be performed. 

