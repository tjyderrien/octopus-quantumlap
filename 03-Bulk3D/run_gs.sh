#!/bin/bash
# Input: $1: SampleLayers, $2: SampleXY
# Purpose: GS calculation

set -e

# ENABLE DRY RUN ? 
DryRun=0
# kPoints=16

source $QuantumLaPsources/machines.env
CheckComputerName
source $QuantumLaPsources/03-Bulk3D/machines.env
echo "[Checking] Variable basePath: $basePath"

# Renaming inputs
SampleLayers=$1
SampleXY=$2
MeshSizeXInput_SI=$3
MeshSizeZInput_SI=$4
kPoints=$5
InitialTemperature=$6
TheoryLevel=$7
Functional=$8

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

  if [ ${Functional} == "mgga_x_tb09" ]; then
   export FunctionalPath="-TB09" 
  fi

  export MeshSizeX_AU=$( python_calc "${MeshSizeX_SI} / ${Bohr2Angstrom}"  ) #Bohr
  export MeshSizeZ_AU=$( python_calc "${MeshSizeZ_SI} / ${Bohr2Angstrom}"  ) #Bohr

  export localPath="${OctopusSlabsRuns}/03-Bulk3D/Si-${SampleXY}x${SampleXY}x${SampleLayers}-k-${kPoints}x${kPoints}x1-Mesh${MeshSizeX_AU}-Mesh${MeshSizeZ_AU}-T${InitialTemperature}${TheoryLevelPath}${FunctionalPath}"
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
  # INP file is modified in root, then copied to the ${localPath} folder  
  echo "[Run_td.sh: ReplaceInputFile] Path for submission: $( pwd ) "
  if [ -e inp ]; then
	  echo ""; echo "Input file INP was renamed to .BAK..."
	  mv inp inp.bak
  fi
  echo "Input file INP was set to ground state."
  cp $OctopusSlabsSources/03-Bulk3D/inp-gs inp || exit

  if [ ${Functional} != "" ]; then
    echo "XCFunctional = ${Functional}" >> inp
  fi
}

BuildSample () {
  # Build the sample
  if [ -e Si.xyz ]; then
    rm Si.xyz
    echo "[Info] Removed the previous sample file."
  fi
  echo "[Checking] AbsPath: $AbsPath"
  echo "[Info] Building the geometry input..."; echo ""
  $OctopusSlabsSources/Samples/SiBulk/build-Si-bulk8.sh >& geometry.in #this is in Angstroms!
  # $AbsPath/Samples/fhi2octopus.sh geometry.in
  echo "[Info] Conversion to octopus format"
  $OctopusSlabsSources/Samples/fhi2octopus.sh geometry.in > Si.xyz #this is now in Bohr!

  # Add the number of lines at the beginning of the file
  numLines=$( cat geometry.in | grep 'atom' | wc -l | awk '{ print $1 }' )

  echo $numLines >> Si.tmp
  echo "# This line is a comment." >> Si.tmp
  cat Si.xyz >> Si.tmp
  cat Si.tmp
  mv Si.tmp Si.xyz

  # Define the missing parameters
  sed -i.bak s/SizeOfTheBox/$SampleXY/g inp
  sed -i.bak s/NumberOfBlocks/$SampleLayers/g inp
#   sed -i.bak s/SizeOfVacuum/$SizeOfVacuum/g inp
  sed -i.bak s/MeshSizeX_AU/$MeshSizeX_AU/g inp
  sed -i.bak s/MeshSizeZ_AU/$MeshSizeZ_AU/g inp
  sed -i.bak s/kPointsValue/${kPoints}/g inp

  if [ ${TheoryLevel} == "independent_particles" ]; then 
        echo "Theory level: independent particles was added to the input file."
        echo "Warning: maybe you meant TDFreezeHxc ? Warning: it is using the same ground state as usual."
        echo "TheoryLevel = independent_particles" >> inp
#        echo "TDFreezeHxc=yes" >> inp #THIS LEAD TO NEVER CONVERGE THE GS. 
  fi

  if [ ${Functional} != "" ]; then
    echo "Theory Functional: you have specified a xc functional. ${Functional} was added to the input file."
    echo "XCFunctional = ${Functional}" >> inp
    echo "FilterPotentials = filter_none" >> inp
  else
    echo "Default LDA functional was used (lda_x + lda_c_pz_mod). See http://octopus-code.org/doc/7.0/html/vars.php?page=alpha." 
  fi
}

SubmitTheJob () {
  echo "== Copying the file in the submission folder... =="
  echo "Source: ${OctopusSlabsRuns}."
  echo "Target: ${localPath}."

  cp ${OctopusSlabsRuns}/03-Bulk3D/inp ${localPath}/
  cp ${OctopusSlabsRuns}/03-Bulk3D/Si.xyz ${localPath}/
  cp $QuantumLaPsources/03-Bulk3D/${SubmitScript} ${localPath}/

  if [ ${DryRun} == 0 ]; then
	  cd ${localPath}
	  ${SubmitCommand} ${SubmitScript} > submit.log
          jobnumber=$( cat submit.log | awk '{ print $4 }')
          cp submit.log $jobnumber.log
          echo "[Submit] Job number: $jobnumber"
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
