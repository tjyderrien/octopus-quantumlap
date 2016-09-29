#!/bin/bash
# Input: $1: SampleLayers, $2: SampleXY
# Purpose: GS calculation

set -e

# ENABLE DRY RUN ? 
DryRun=1
kPoints=16

source machines.env

# Renaming inputs
SampleLayers=$1
SampleXY=$2
# VacuumSizeInput=$3
MeshSizeXInput=$4
MeshSizeZInput=$5

python_calc () {
    python -c "print ($@)"
}

# Conversion factors
# Bohr2Angstrom=0.52917721067 #AngstromPerBohr

# Renaming inputs
ProcessArguments () { # Conversion factors
  Bohr2Angstrom=0.52917721067 #AngstromPerBohr

  # Renaming inputs
  SampleLayers=$1
  SampleXY=$2
  # SizeOfVacuumAng=${VacuumSizeInput} #10. #Angstrom
  # SizeOfVacuum=$( python_calc "${SizeOfVacuumAng} / ${Bohr2Angstrom}" ) #Bohr
  # SizeOfVacuumName=$( python_calc "int(round( ${SizeOfVacuumAng} ))" )
  MeshSizeXAng=${MeshSizeXInput} #0.23 #Angstrom
  MeshSizeZAng=${MeshSizeZInput} #should be smaller, as we have a surface
  export MeshSizeX=$( python_calc "${MeshSizeXAng} / ${Bohr2Angstrom}"  ) #Bohr
  export MeshSizeZ=$( python_calc "${MeshSizeZAng} / ${Bohr2Angstrom}"  ) #Bohr

  export localPath="${basePath}/Si-${SampleXY}x${SampleXY}x${SampleLayers}-k-${kPoints}x${kPoints}x1-Mesh${MeshSizeX}-Mesh${MeshSizeZ}"
  echo "== BUILDING FOLDER ... ${localPath} =="
  mkdir -p ${localPath}

}

print_usage () {
  if [ -e ${SampleXY} ]; then
	  echo "Usage: "
	  echo "./run.sh <Number of Z blocks> <Number of XY atoms> <Vacuum Size> <MeshXY size, MeshZ size>"
	  exit
  else
	  echo "** Sample size: Z=$SampleLayers layers. XY=${SampleXY} atoms.**"
          echo "** Mesh size: Z=${MeshSizeZInput}. X=${MeshSizeXInput}. **"
  fi
}

ReplaceInputFile () {
  # We would like to start the ground state
  if [ -e ${localPath}/inp ]; then
	  echo "Input file INP was renamed to .BAK..."
	  mv ${localPath}/inp ${localPath}/inp.bak
  fi
  echo "Input file INP was set to ground state."
  cp ${basePath}/inp-gs ${localPath}/inp || exit
}

BuildSample () {
  # Build the sample
  if [ -e ${localPath}/Si.xyz ]; then
    rm ${localPath}/Si.xyz
  fi
  $AbsPath/Samples/SiBulk/build-Si-bulk8.sh >& geometry.in
  # $AbsPath/Samples/fhi2octopus.sh geometry.in 
  $AbsPath/Samples/fhi2octopus.sh geometry.in > ${localPath}/Si.xyz

  # Add the number of lines at the beginning of the file
  numLines=$( cat geometry.in | grep 'atom' | wc -l | awk '{ print $1 }' )

  echo $numLines >> ${localPath}/Si.tmp
  echo "# This line is a comment." >> ${localPath}/Si.tmp
  cat ${localPath}/Si.xyz >> ${localPath}/Si.tmp
  cat ${localPath}/Si.tmp
  mv ${localPath}/Si.tmp ${localPath}/Si.xyz

  # Define the missing parameters
  sed -i.bak s/SizeOfTheBox/$SampleXY/g ${localPath}/inp
  sed -i.bak s/NumberOfBlocks/$SampleLayers/g ${localPath}/inp
  sed -i.bak s/SizeOfVacuum/$SizeOfVacuum/g ${localPath}/inp
  sed -i.bak s/MeshSizeX/$MeshSizeX/g ${localPath}/inp
  sed -i.bak s/MeshSizeZ/$MeshSizeZ/g ${localPath}/inp
  sed -i.bak s/kPointsValue/${kPoints}/g ${localPath}/inp
}

SubmitTheJob () {
  cp ${basePath}/${SubmitScript} ${localPath}/
  if [ ${DryRun} == 0 ]; then
	  cd ${localPath}
	  ${SubmitCommand} ${SubmitScript}
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
ReplaceInputFile 
BuildSample
SubmitTheJob
