#!/bin/bash
# Input: $1: SampleLayers, $2: SampleXY
# Purpose: GS calculation

set -e

# ENABLE DRY RUN ? 
DryRun=1
kPoints=1

source $QuantumLaPsources/machines.env
CheckComputerName
source $QuantumLaPsources/02-Slab2D/machines.env
echo "[Checking] Variable basePath: $basePath"

# Renaming inputs (Input units: Angstroms)
SampleLayers=$1
SampleXY=$2
VacuumSizeInput_SI=$3
MeshSizeXInput_SI=$4
MeshSizeZInput_SI=$5

TheoryLevel=$6
Functional=$7

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
  SizeOfVacuum_SI=${VacuumSizeInput_SI} #10. #Angstrom
  SizeOfVacuum_AU=$( python_calc "${SizeOfVacuum_SI} / ${Bohr2Angstrom}" ) #Bohr
  SizeOfVacuumName_SI=$( python_calc "int(round( ${SizeOfVacuum_SI} ))" )
  MeshSizeX_SI=${MeshSizeXInput_SI} #0.23 #Angstrom
  MeshSizeZ_SI=${MeshSizeZInput_SI} #should be smaller, as we have a surface
  
  if [ ${TheoryLevel} == "independent_particles" ]; then
    export TheoryLevelPath="-Indep_particles"
  else
    export TheoryLevelPath=""
  fi

  case "${Functional}" in
    "TB09")
      FunctionalPath="-TB09"
      ;;
    "LDA")
      FunctionalPath=""
      ;;
    *)
      echo "This functional was not planned. Modify run_td.sh."
      exit 1
      ;;
  esac
  
  export MeshSizeX_AU=$( python_calc "${MeshSizeX_SI} / ${Bohr2Angstrom}"  ) #Bohr
  export MeshSizeZ_AU=$( python_calc "${MeshSizeZ_SI} / ${Bohr2Angstrom}"  ) #Bohr

  export localPath="${OctopusSlabsRuns}/02-Slab2D/Si-${SampleXY}x${SampleXY}x${SampleLayers}-k-${kPoints}x${kPoints}x1-L${SizeOfVacuumName_SI}-Mesh${MeshSizeX_AU}-Mesh${MeshSizeZ_AU}-${TheoryLevelPath}${FunctionalPath}" #Mesh: stored in Bohr in the path. #Vacuum: stored in angstroms. 
  echo "== BUILDING FOLDER ... ${localPath} =="
  mkdir -p ${localPath}

}

print_usage () {
  if [ -e ${SampleXY} ]; then
	  echo "Usage: "
	  echo "./run.sh <Number of Z blocks> <Number of XY atoms> <Vacuum Size (SI, angstroms) > <MeshXY size (SI, angstroms), MeshZ size (SI, angstroms)>"
	  exit
  else
	  echo "** Sample size: Z=$SampleLayers layers. XY=${SampleXY} atoms.**"
          echo "** Mesh size (SI, angstroms): Z=${MeshSizeZInput_SI}. X=${MeshSizeXInput_SI}. **"
  fi
}

ReplaceInputFile () {
  # We would like to start the ground state
  if [ -e inp ]; then
	  echo "Input file INP was renamed to .BAK..."
	  mv inp inp.bak
  fi
  echo "Input file INP was set to ground state."
  cp $OctopusSlabsSources/02-Slab2D/inp-gs inp || exit
  
  if [ ${Functional} == "TB09" ]; then
    echo "XCFunctional = mgga_x_tb09 + lda_c_pw" >> inp
    echo "FilterPotentials = filter_none" >> inp
  fi

  if [ ${Functional} == "LDA" ]; then
    echo "XCFunctional = lda_x + lda_c_pz_mod" >> inp
  fi
}

BuildSample () {
  # Build the sample [Si.xyz is in SI (angstroms)]. 
  if [ -e Si.xyz ]; then
    rm Si.xyz
    echo "[Info] Removed the previous sample file."
  fi
  echo "[Checking] AbsPath: $AbsPath"
  echo "[Info] Building the geometry input..."; echo ""
  $OctopusSlabsSources/Samples/Si100-H/write-geom.py ${SampleLayers} ${SampleXY} >& geometry.in
  $OctopusSlabsSources/Samples/fhi2octopus.sh geometry.in > Si.xyz #Conversion to Bohrs

  # Add the number of lines at the beginning of the file
#   numLines=$( wc -l ${localPath}/Si.xyz | awk '{ print $1 }' )
  numLines=$( cat geometry.in | grep 'atom' | wc -l | awk '{ print $1 }' )

  echo $numLines >> Si.tmp
  echo "# This line is a comment." >> Si.tmp
  cat Si.xyz >> Si.tmp
  cat Si.tmp
  mv Si.tmp Si.xyz

  # Define the missing parameters in input file (always in atomic units)
  sed -i.bak s/SizeOfTheBox/$SampleXY/g inp #number of layers
  sed -i.bak s/NumberOfBlocks/$SampleLayers/g inp
  sed -i.bak s/SizeOfVacuum_AU/$SizeOfVacuum_AU/g inp
  sed -i.bak s/MeshSizeX_AU/$MeshSizeX_AU/g inp
  sed -i.bak s/MeshSizeZ_AU/$MeshSizeZ_AU/g inp
  sed -i.bak s/kPointsValue/${kPoints}/g inp
  
  if [ ${TheoryLevel} == "independent_particles" ]; then 
      echo "Theory level: independent particles was added to the input file."
      echo "Warning: maybe you meant TDFreezeHxc ? Warning: it is using the same ground state as usual."
      echo "TheoryLevel = independent_particles" >> inp
#      echo "TDFreezeHxc=yes" >> inp #THIS LEAD TO NEVER CONVERGE THE GS. 
  fi
}

SubmitTheJob () {
  echo "== Copying the file in the submission folder... =="
  echo "Source: ${OctopusSlabsRuns}."
  echo "Target: ${localPath}."
  
  cp ${OctopusSlabsRuns}/02-Slab2D/inp ${localPath}/
  cp ${OctopusSlabsRuns}/02-Slab2D/Si.xyz ${localPath}/
  cp ${QuantumLaPsources}/02-Slab2D/${SubmitScript} ${localPath}/
  
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
