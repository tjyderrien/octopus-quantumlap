#!/bin/bash
# Input: $1: SampleLayers, $2: SampleXY
# Purpose: GS calculation

set -e

# ENABLE DRY RUN ? 
DryRun=0
# kPoints=16

source $QuantumLaPsources/octopus_slabs/machines.env
CheckComputerName
LoadEnvLibraries
source $QuantumLaPsources/octopus_slabs/machines.env
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
Mission=$9
NumNodes=${10}
AngleIncidence=${11}


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


#  export localPath="${QuantumLaPruns}/03-Bulk3D/Si-${SampleXY}x${SampleXY}x${SampleLayers}-k-${kPoints}x${kPoints}x1-Mesh${MeshSizeX_AU}-Mesh${MeshSizeZ_AU}-T${InitialTemperature}${TheoryLevelPath}${FunctionalPath}"

  export localPath="${QuantumLaPruns}/03-Bulk3D/Si-${SampleXY}x${SampleXY}x${SampleLayers}-k-${kPoints}x${kPoints}x1-Angle${AngleIncidence}-Mesh${MeshSizeX_AU}-Mesh${MeshSizeZ_AU}-T${InitialTemperature}${TheoryLevelPath}${FunctionalPath}"

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
  cp $QuantumLaPsources/03-Bulk3D/inp-gs inp || exit

  if [ ${Functional} == "TB09" ]; then
    echo "XCFunctional = mgga_x_tb09 + lda_c_pw" >> inp
    echo "FilterPotentials = filter_none" >> inp
  fi

  if [ ${Functional} == "LDA" ]; then
    echo "XCFunctional = lda_x + lda_c_pz_mod" >> inp
  fi


}

BuildSample () {
  # Build the sample
  if [ -e Si.xyz ]; then
    rm Si.xyz
    echo "[Info] Removed the previous sample file."
  fi
  echo "[Checking] AbsPath: $AbsPath"
  # echo "[Info] Building the geometry input..."; echo ""
  # $QuantumLaPsources/Samples/SiBulk/build-Si-bulk8.sh >& geometry.in #this is in Angstroms!
  # $AbsPath/Samples/fhi2octopus.sh geometry.in
  # echo "[Info] Conversion to octopus format"
  # $QuantumLaPsources/Samples/fhi2octopus.sh geometry.in > Si.xyz #this is now in Bohr!

  # Add the number of lines at the beginning of the file
  # numLines=$( cat geometry.in | grep 'atom' | wc -l | awk '{ print $1 }' )

  # echo $numLines >> Si.tmp
  # echo "# This line is a comment." >> Si.tmp
  # cat Si.xyz >> Si.tmp
  # cat Si.tmp
  # mv Si.tmp Si.xyz

  # Define the missing parameters
  sed -i.bak s/SizeOfTheBox/$SampleXY/g inp
  sed -i.bak s/NumberOfBlocks/$SampleLayers/g inp
#   sed -i.bak s/SizeOfVacuum/$SizeOfVacuum/g inp
  sed -i.bak s/MeshSizeX_AU/$MeshSizeX_AU/g inp
  sed -i.bak s/MeshSizeZ_AU/$MeshSizeZ_AU/g inp
  sed -i.bak s/kPointsValue/${kPoints}/g inp
  sed -i.bak s/theta1_inp/${AngleIncidence}/g inp

  if [ ${TheoryLevel} == "independent_particles" ]; then 
        echo "Theory level: independent particles was added to the input file."
        echo "Warning: maybe you meant TDFreezeHxc ? Warning: it is using the same ground state as usual."
        echo "TheoryLevel = independent_particles" >> inp
#        echo "TDFreezeHxc=yes" >> inp #THIS LEAD TO NEVER CONVERGE THE GS. 
  fi

}

SubmitTheJob () {
  echo "== Copying the file in the submission folder... =="
  echo "Source: ${QuantumLaPruns}."
  echo "Target: ${localPath}."

  cp -u inp ${localPath}/ #just update if file was not found, or is older.
  # cp ${QuantumLaPsources}/03-Bulk3D/Si.xyz ${localPath}/ #no longer needed, as it is included in input file
  cp ${QuantumLaPsources}/03-Bulk3D/${SubmitScript} ${localPath}/

  sed -i.bak s/Mission/${Mission}/g ${localPath}/${SubmitScript} 
## TODO: this region is machine dependent
  sed -i.bak s/NumNodes/${NumNodes}/g ${localPath}/${SubmitScript} 
## 

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
