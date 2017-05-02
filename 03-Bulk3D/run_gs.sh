#!/bin/bash
# Input: $1: SampleLayers, $2: SampleXY
# Purpose: GS calculation

set -e

# ENABLE DRY RUN ? 
DryRun=0
# kPoints=16

source machines.env

# Renaming inputs
SampleLayers=$1
SampleXY=$2
MeshSizeXInput_SI=$3
MeshSizeZInput_SI=$4
kPoints=$5
InitialTemperature=$6
TheoryLevel=$7

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
  MeshSizeX_SI=${MeshSizeXInput_SI} #0.23 #Angstrom
  MeshSizeZ_SI=${MeshSizeZInput_SI} #should be smaller, as we have a surface
  if [ ${TheoryLevel} == "independent_particles" ]; then
    export TheoryLevelPath="-Indep_particles"
  else
    export TheoryLevelPath=""
  fi

  export MeshSizeX_AU=$( python_calc "${MeshSizeX_SI} / ${Bohr2Angstrom}"  ) #Bohr
  export MeshSizeZ_AU=$( python_calc "${MeshSizeZ_SI} / ${Bohr2Angstrom}"  ) #Bohr

  export localPath="${basePath}/Si-${SampleXY}x${SampleXY}x${SampleLayers}-k-${kPoints}x${kPoints}x1-Mesh${MeshSizeX_AU}-Mesh${MeshSizeZ_AU}-T${InitialTemperature}${TheoryLevelPath}"
  echo "== BUILDING FOLDER ... ${localPath} =="
  mkdir -p ${localPath}

}

print_usage () {
  if [ -e ${SampleXY} ]; then
	  echo "Usage: "
	  echo "./run.sh <Number of Z blocks> <Number of XY atoms> <MeshXY size, MeshZ size>"
	  exit
  else
	  echo "** Sample size: Z=$SampleLayers layers. XY=${SampleXY} atoms.**"
          echo "** Mesh size: Z=${MeshSizeZInput_SI}. X=${MeshSizeXInput_SI}. **"
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
  $AbsPath/Samples/SiBulk/build-Si-bulk8.sh >& geometry.in #this is in Angstroms!
  # $AbsPath/Samples/fhi2octopus.sh geometry.in 
  $AbsPath/Samples/fhi2octopus.sh geometry.in > ${localPath}/Si.xyz #this is now in Bohr!

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
#   sed -i.bak s/SizeOfVacuum/$SizeOfVacuum/g ${localPath}/inp
  sed -i.bak s/MeshSizeX_AU/$MeshSizeX_AU/g ${localPath}/inp
  sed -i.bak s/MeshSizeZ_AU/$MeshSizeZ_AU/g ${localPath}/inp
  sed -i.bak s/kPointsValue/${kPoints}/g ${localPath}/inp

  if [ ${TheoryLevel} == "independent_particles" ]; then
        echo "Theory level: independent particles was added to the input file."
        echo "TheoryLevel = independent_particles" >> inp
  fi
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
