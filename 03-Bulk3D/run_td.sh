#!/bin/bash
# Input: $1: SampleLayers, $2: SampleXY
# Purpose: TD-evolution. 

# set -e
source machines.env
source functions.sh
# ENABLE DRY RUN ? 
DryRun=1
MoveIons=0 #cost: 3x runtime(MoveIons=no)
kPoints=8 # 8: convergence of density is fine, 24: convergence of currents is fine

# Renaming inputs to replace into the INP file.
SampleLayers=$1
SampleXY=$2
MeshSizeXInput_SI=$3
MeshSizeZInput_SI=$4
lambda1=$5
lambda2=$6
Efield1=$7
Efield2=$8
tau1=$9
tau2=${10}
CEP1=${11}
CEP2=${12}
delay=${13}
InitialTemperature=${14}
PulseShape=${15}
InducedField=${16}

ProcessArguments () { # Conversion factors for sample
  Bohr2Angstrom=0.52917721067 #AngstromPerBohr

  # Renaming inputs
  SampleLayers=$1
  SampleXY=$2
#   SizeOfVacuum_SI=${VacuumSizeInput} #Angstrom
#   SizeOfVacuum_AU=$( python_calc "${SizeOfVacuum_SI} / ${Bohr2Angstrom}" ) #Bohr
#   SizeOfVacuumName_SI=$( python_calc "int(round( ${SizeOfVacuum_SI} ))" )
  #MeshSizeAng=0.23 #Angstrom
  export MeshSizeX_SI=${MeshSizeXInput_SI} #0.23 #Angstrom
  export MeshSizeZ_SI=${MeshSizeZInput_SI} #should be smaller, as we have a surface
  export MeshSizeX_AU=$( python_calc "${MeshSizeX_SI} / ${Bohr2Angstrom}"  ) #Bohr
  export MeshSizeZ_AU=$( python_calc "${MeshSizeZ_SI} / ${Bohr2Angstrom}"  ) #Bohr
#  MeshSize=$( python_calc "${MeshSizeAng} / ${Bohr2Angstrom}"  ) #Bohr
}

DefineParameters () {
  # Input laser parameters
  if [ $MoveIons == 1 ]; then
    RelaxationTime=$( python_calc "500e-15"  )
  else
    RelaxationTime=$( python_calc "60e-15" )
  fi

# echo "Generating pulses..."

#  echo "TODO: for now, the Python laser generator only makes nice pictures."
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
  echo "[Run_TD] Pulses generated."
  
  if [ $InducedField == 1 ]; then
     GaugeFieldValue="yes"
     GaugeFieldText="-InducedFields"
  else
     GaugeFieldValue="no"
     GaugeFieldText=""
  fi

  MoveIonsChain=""
  if [ ${MoveIons} == 1 ]; then
	MoveIonsChain="-MoveIons"
  else
        MoveIonsChain=""
  fi
  export groundStatePath=${basePath}/Si-${SampleXY}x${SampleXY}x${SampleLayers}-k-${kPoints}x${kPoints}x1-Mesh${MeshSizeX_AU}-Mesh${MeshSizeZ_AU}-T${InitialTemperature}
  export localPath=${groundStatePath}/"lambda"${lambda1}-"tau"${tau1}-"field"${Efield1}-"CEP"${CEP1}-${PulseShape}/"lambda"${lambda2}-"tau"${tau2}-"field"${Efield2}-"CEP"${CEP2}-"delay"${delay}-${PulseShape}${GaugeFieldText}${MoveIonsChain}
  if [ ! -d ${localPath} ]; then
	mkdir -p ${localPath}
  fi
}

print_usage () {
  if [ -e ${SampleXY} ]; then
    echo "Usage: "
    echo "./run.sh <Number of Z blocks> <Number of XY atoms> <Mesh size XY (Angstrom), Mesh size Z (Ang.)> <lambda1, lambda2 (m)>, <Efield1, Efield2 (V/m)>, <tau1, tau2 (s)>, <CEP1 (rad), CEP2 (rad)>, <delay (s)>"
    exit
  else
    echo ""
    echo "[Info] ** Sample size: Z=${SampleLayers} layers. XY=${SampleXY} atoms.**"
    echo "[Info] ** Mesh size: Z=${MeshSizeZInput_SI}. X=${MeshSizeXInput_SI}. **"
    echo ""
  fi
}

ReplaceInputFile () {
  # Define the material parameters
  # We would like to start the td-calculation
  # INP file is modified in root, then copied to the ${localPath} folder  
  if [ -e inp ]; then
	  echo ""; echo "Input file INP was renamed to .BAK..."
	  mv inp inp.bak
  fi
  echo "Input file INP was set to TIME-DEPENDENT RUN"
  cp inp-td inp
  echo ""; echo "REPLACING MATERIAL PARAMETERS..."
  sed -i.bak s/SizeOfTheBox/${SampleXY}/g inp
  sed -i.bak s/NumberOfBlocks/${SampleLayers}/g inp
#   sed -i.bak s/SizeOfVacuum_AU/${SizeOfVacuum_AU}/g inp
  sed -i.bak s/MeshSizeX_AU/${MeshSizeX_AU}/g inp
  sed -i.bak s/MeshSizeZ_AU/${MeshSizeZ_AU}/g inp
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
  sed -i.bak s/FieldShapeA/FieldColorA-${PulseShape}/g inp
  sed -i.bak s/FieldShapeB/FieldColorB-${PulseShape}/g inp
  sed -i.bak s/ReplaceMeGaugeField/${GaugeFieldValue}/g inp

  if [ ${MoveIons} == 1 ]; then
	echo "MoveIons? Addind the function in Input file..."
	echo "MoveIons = yes" >> inp
  fi

}

BuildSample() {
  # Build the sample
  if [ -e Si.xyz ]; then
    rm Si.xyz
    echo "[Info] Removed the previous sample file."
  fi
  echo "[Info] Building the geometry input..."; echo ""
  $AbsPath/Samples/SiBulk/build-Si-bulk8.sh >& geometry.in #build the geometry.in file.
  echo "[Info] Conversion to octopus format" 
  $AbsPath/Samples/fhi2octopus.sh geometry.in >& Si.xyz #Convert to octopus format. 

  # Add the number of lines at the beginning of the file
  numLines=$( wc -l Si.xyz | awk '{ print $1 }' )

  echo $numLines >> Si.tmp
  echo "# This line is a comment." >> Si.tmp
  cat Si.xyz >> Si.tmp
  echo ""; echo "** Content of the Si.xyz file."
  cat Si.tmp
  mv Si.tmp Si.xyz
}

SubmitTheJob () {
  echo "== CHECKING GROUND STATE INTEGRITY =="
  #If ground state exists, we shall copy the ground state in the target folder!
  #else: don't run the ground, just stop the program and tell it. 

  if [ $DryRun == 1 ]; then
	DryRunChain="--dry-run"
  else
	DryRunChain=""
  fi
  if [ -e ${groundStatePath}-Save ]; then
          echo "[Copy] Source: ${groundStatePath}-Save"
          echo "[Copy] Target: ${groundStatePath}"
          ./DuplicateGroundState.sh ${groundStatePath}-Save/ ${groundStatePath}/ $DryRunChain
  else
    echo "GROUND STATE DOES NOT EXIST."
    echo "Expected path: ${groundStatePath}-Save"
    echo "Run it manually using run_gs.sh ${SampleLayers} ${SampleXY}."
    exit
  fi 
  echo ""
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

