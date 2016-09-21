#!/bin/bash

echo "=== LIST OF AVAILABLE JOBS ==="
source machines.env
CheckComputerName
source WhoAreTheRunningJobs.sh

echo "==== CHECK FOR DUPLICATES ==="
./WhoAreTheRunningJobs.sh | uniq -c 

echo "==== CALCULATION OF EFFICIENCIES ===="

for jobID in $( GetJobID )
do
	#Get the corresponding log file
	logFile=$( find | grep "$jobID" | grep 'log' | grep 'octopus' )
	echo "LogFile=$logFile"
	# Catch the latest lines with TD data. 
	#### COMPUTE THE PHYSICAL TIME ####
	ElapsedPhysicalTime=$( tail -n1 $logFile | awk '{ print $2*1E15*2.44e-17 }' ) #s
	#HumanTimePerLoop=$( tail -n1 $logFile | awk '{ print $5 }' ) #s

	### COMPUTE THE ELAPSED HUMAN TIME ###
	# Method1: uses squeue, but limited to 24 hours. 
	ExtractTheElapsedTime=$( $QueueCommand | grep "$jobID" | awk "{ print ${RunningTimeColumn} }" | tr ':' ' ' | awk '{ print $1 }' )
#	echo "ExtractedTime = $ExtractTheElapsedTime"
	if [ -z ${ExtractTheElapsedTime} ]; then #TODO: find out why this test.
	#2 digits as < 1 h yet
		ElapsedHumanTime=$( $QueueCommand | grep "$jobID" | awk '{ print $6 }' | tr ':' ' ' | awk '{ print ($2+$1*60)/(60.*60.) }' ) #in hours
	else 
		ElapsedHumanTime=$( $QueueCommand | grep "$jobID" | awk '{ print $6 }' | tr ':' ' ' | awk '{ print ($3+$2*60.+$1*3600.)/(60.*60.) }' ) #in hours
	fi

#TODO: Method2: uses ls -lhtr octopus-td*.log"

	#### COMPUTE THE EFFICIENCY ###
	Efficiency=$( python -c "print $ElapsedPhysicalTime/$ElapsedHumanTime * 24." )

	echo "Elapsed Physical Time: ${ElapsedPhysicalTime} fs"
	echo "Elapsed Human Time: ${ElapsedHumanTime} hours"
	echo "=== EFFICIENCY: $Efficiency fs/day. ==="
done


#
#squeue | grep 'tjyderri'


