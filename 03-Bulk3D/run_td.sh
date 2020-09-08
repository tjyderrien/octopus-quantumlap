#!/bin/bash
# Input: $1: SampleLayers, $2: SampleXY
# Purpose: TD-evolution. 

set -e
source $QuantumLaPsources/machines.env
CheckComputerName
source $QuantumLaPsources/03-Bulk3D/machines.env
echo "[Checking] Variable basePath: $basePath"

source $OctopusSlabsSources/03-Bulk3D/functions.sh
# ENABLE DRY RUN ? 
DryRun=0
MoveIons=0 #cost: 3x runtime(MoveIons=no), and ForceZeroForce is not enabled!
kPoints=24 # 8: convergence of density is fine, 24: convergence of currents is fine

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
CrossPolarization=${17}
TheoryLevel=${18}
Functional=${19}
FreezeXC=${20}
RiseCycle=${21}
Mission=${22}
NumNodes=${23}
GaugeVectorFieldAmplitude=${24}
GaugeFieldDelayValue=${25}
UseLinearResponse=${26}

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
    export RelaxationTime=$( python_calc "500e-15"  )
  else
    export RelaxationTime=$( python_calc "45e-15" )
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
#  ComputingTime=$( python_calc " float( (8.*${tau1}/${sigmaTau} + ${delay} + ${RelaxationTime}) ) " )
  ComputingTime=$( python_calc " float ( ( 8./5. * ( 4*${lambda1}/299792458.0 +  ${tau1}) ) ) " )
#   tauRiseA=$( python_calc " float( ${RiseCycle} / (2.* 3.141592653 * 2.99792458e8 / ${lambda1} ) ) "  )
#   tauRiseB=$( python_calc " float( ${RiseCycle} / (2.* 3.141592653 * 2.99792458e8 / ${lambda2} ) ) "  )
  if [ ${PulseShape} == "StH4" ]; then
    SinPower=4
    echo "** Info: power of Soft Top Hat enveloppe is defined to 4."
  else 
    SinPower=2
    echo "** Info: power of Soft Top Hat enveloppe is defined to 2."
  fi
# CONCEPT OF THE PULSE CAPTURE
# NOTE: flood_td.sh generates the list of pulses we need. 
# NOTE: run_td.sh addresses each of the requests to octopus/Pulses/python script who calls spp-extended-theory/Keldysh/libKeldyshPulses.py:PulseSquaredSinTemporalShape() to generates the envelope. 
# NOTE: The envelope is saved into a 2-columns file that will be used by octopus tdf_from_file. 
# NOTE: Units are in atomic units. 
# NOTE: run_td.sh puts the result in a file named PulseA|PulseB.dat in the execution folder. 
  $QuantumLaPsources/Pulses/GeneratePulse.py ${tau1} ${RiseCycle} ${Efield1} ${lambda1} ${CEP1} 0. ${ComputingTime} ${SinPower}
  mv Pulses.dat PulseA.dat
  $QuantumLaPsources/Pulses/GeneratePulse.py ${tau2} ${RiseCycle} ${Efield2} ${lambda2} ${CEP2} ${delay} ${ComputingTime} ${SinPower}
  mv Pulses.dat PulseB.dat
  
  echo "[Run_TD] Pulses generated."
  
  if [ ${CrossPolarization} == 1 ]; then
     CrossPolarizationText="-CrossPola"
     CrossPolaAngle="0.5*pi"
  else
     CrossPolaAngle="0."
  fi

  if [ $InducedField == 1 ]; then
     GaugeFieldValue="yes"
     GaugeFieldText="-InducedFields"
  else
     GaugeFieldValue="no"
     GaugeFieldText=""
  fi

# PROBE INPUT
#  if [ ${GaugeVectorFieldAmplitude} == 0 ]; then #BUG: RISKY

  if [ ${UseLinearResponse} == 0 ]; then
    GaugeVectorFieldText=""
  else
    GaugeVectorFieldText="-Probe${GaugeVectorFieldAmplitude}-ProbeDelay${GaugeFieldDelayValue}"
  fi

  MoveIonsChain=""
  if [ ${MoveIons} == 1 ]; then
	MoveIonsChain="-MoveIons"
  else
        MoveIonsChain=""
  fi

  if [ ${TheoryLevel} == "independent_particles" ]; then
    TheoryLevelPath="-Indep_particles"
  else
    TheoryLevelPath=""
  fi

  if [ ${FreezeXC} == 1 ]; then 
    FreezeXCchain="-FreezeXC"
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

  echo "Generating the folder name..."
  export groundStatePath=${OctopusSlabsRuns}/03-Bulk3D/Si-${SampleXY}x${SampleXY}x${SampleLayers}-k-${kPoints}x${kPoints}x1-Mesh${MeshSizeX_AU}-Mesh${MeshSizeZ_AU}-T${InitialTemperature}${TheoryLevelPath}${FunctionalPath}
  export localPath=${groundStatePath}/"lambda"${lambda1}-"tau"${tau1}-"field"${Efield1}-"CEP"${CEP1}-${PulseShape}-SinPower${SinPower}-${RiseCycle}${FreezeXCchain}${GaugeVectorFieldText}/"lambda"${lambda2}-"tau"${tau2}-"field"${Efield2}-"CEP"${CEP2}-"delay"${delay}-${PulseShape}${CrossPolarizationText}${GaugeFieldText}${MoveIonsChain}
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
  # INP file is modified in the octopus-slabs root, then copied to the ${localPath} folder  
  echo "[Run_td.sh: ReplaceInputFile] Path for submission: $( pwd ) "
  if [ -e inp ]; then
	  echo ""; echo "Input file INP was renamed to .BAK..."
	  mv inp inp.bak
  fi
  echo "Input file INP was set to TIME-DEPENDENT RUN"
  cp $OctopusSlabsSources/03-Bulk3D/inp-td inp
  echo ""; echo "REPLACING MATERIAL PARAMETERS for $localPath"
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
  sed -i.bak s/ReplaceMeTheta/${CrossPolaAngle}/g inp
  if [ ${UseLinearResponse} == 1 ]; then
echo "# Linear response probing" >> inp
echo "%GaugeVectorField #useful for probing" >> inp
echo "   GaugeVectorFieldAmplitude | 0 | 0  #to be replaced by Script" >> inp
echo "%" >> inp
echo "GaugeFieldDelaySI = GaugeFieldDelayValue" >> inp
echo "GaugeFieldDelay  = GaugeFieldDelaySI / TimeAUtoSI" >> inp
    sed -i.bak s/GaugeVectorFieldAmplitude/${GaugeVectorFieldAmplitude}/g inp #THIS GENERATES HUGE BUG if used for usual TD propagation
    sed -i.bak s/GaugeFieldDelayValue/${GaugeFieldDelayValue}/g inp

  fi

  if [ ${MoveIons} == 1 ]; then
	echo "MoveIons? Addind the function in Input file..."
	echo "MoveIons = yes" >> inp
  fi

  if [ ${TheoryLevel} == "independent_particles" ]; then
        echo "Theory level: independent particles was added to the input file."
        echo "WARNING: don't mix with TDFreezeXC !"
        echo "TheoryLevel = independent_particles" >> inp 
  fi

  if [ ${FreezeXC} == 1 ]; then
    echo "TDFreezeHXC = yes" >> inp    
  fi

  if [ ${Functional} == "TB09" ]; then
    echo "XCFunctional = mgga_x_tb09 + lda_c_pw" >> inp
    sed -i 's/TDPropagator/#TDPropagator/g' inp # disable aetrs
    echo "TDPropagator = etrs" >> inp #and replace by etrs only. 
    echo "TDStepsWithSelfConsistency = all_steps" >> inp #with self-consistent
    echo "FilterPotentials = filter_none" >> inp
  fi

  if [ ${Functional} == "LDA" ]; then
    echo "XCFunctional = lda_x + lda_c_pz_mod" >> inp 
  fi
  # cd - >& /dev/null

}

BuildSample() {
  # Build the sample
  if [ -e Si.xyz ]; then
    rm Si.xyz
    echo "[Info] Removed the previous sample file."
  fi
  echo "[Checking] AbsPath: $AbsPath"
  echo "[Info] Building the geometry input..."; echo ""
  $OctopusSlabsSources/Samples/SiBulk/build-Si-bulk8.sh >& geometry.in #build the geometry.in file.
  echo "[Info] Conversion to octopus format" 
  $OctopusSlabsSources/Samples/fhi2octopus.sh geometry.in >& Si.xyz #Convert to octopus format. 

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

# Checking if save of ground state exists. #NOTE: this made sense when it was long and hard to compute ground states. 
#  if [ -e ${groundStatePath}-Save ]; then
#          echo "[Copy] Source: ${groundStatePath}-Save"
#          echo "[Copy] Target: ${groundStatePath}"
#          # $OctopusSlabsSources/02-Slab2D/DuplicateGroundState.sh ${groundStatePath}-Save/ ${groundStatePath}/ $DryRunChain
#  else
#    echo "SAVE OF GROUND STATE DOES NOT EXIST."
#    echo "Expected path: ${groundStatePath}-Save"
#    exit
#  fi

 
  echo ""
  echo "== Copying the file in the submission folder... =="
  echo "Source: ${OctopusSlabsRuns}."
  echo "Target: ${localPath}."
#  cp ${OctopusSlabsRuns}/03-Bulk3D/inp ${localPath}/ #induces an error: parameters were not replaced systematically
# Good path is actually $( pwd ).
  cp inp ${localPath}/ && echo "Octopus INP source successfully copied to destination."
  cp PulseA.dat ${localPath}/ && echo "Successfully copied Pulse A to destination."
  cp PulseB.dat ${localPath}/ && echo "Successfully copied Pulse B to destination."
#  cp ${basePath}/submit_mpi.sge ${localPath}/
  cp ${OctopusSlabsRuns}/03-Bulk3D/Si.xyz ${localPath}/
  
  # Update the Submission script with the mission name
  cp $QuantumLaPsources/03-Bulk3D/$SubmitScript ${localPath}/ && echo "Successfully copied SLURM script to destination."
  sed -i.bak s/Mission/${Mission}/g ${localPath}/${SubmitScript}
  sed -i.bak s/NumNodes/${NumNodes}/g ${localPath}/${SubmitScript}

  # Start the run
  if [ ${DryRun} == 0 ]; then
          cd ${localPath}
          ${SubmitCommand} ${SubmitScript} > ${localPath}/submit.log
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
DefineParameters
ReplaceInputFile #this is performed in the common directory!
BuildSample #also performed in directory
# ReplaceInputFile 
SubmitTheJob #then the file are move in the directory where calculation will be performed. 

