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
	logDir=$( dirname  $( find | grep "$jobID" | tail -n1 ) )
	logFile="$logDir/$ResultFile"
	echo "LogFile=$logFile"
	# Catch the latest lines with TD data. 
	#### COMPUTE THE PHYSICAL TIME ####
	ElapsedPhysicalTime=$( grep 'SCF CYCLE' $logFile | tail -n1 | awk '{ print $6 }' ) #Iteration number
	#HumanTimePerLoop=$( tail -n1 $logFile | awk '{ print $5 }' ) #s

	### COMPUTE THE ELAPSED HUMAN TIME ###
	# Method1: uses squeue, but limited to 24 hours. 
	ExtractTheElapsedTime=$( $QueueCommand | grep "$jobID" | awk "{ print ${RunningTimeColumn} }" | tr ':' ' ' | awk '{ print $1 }' )
#	echo "ExtractedTime = $ExtractTheElapsedTime"
	if [ -z ${ExtractTheElapsedTime} ]; then #TODO: change the method using date -u -d, and compare durations with 1st Jan 2000. 
	#2 digits as < 1 h yet
		ElapsedHumanTime=$( $QueueCommand | grep "$jobID" | awk '{ print $6 }' | tr ':' ' ' | awk '{ print ($2+$1*60)/(60.*60.) }' ) #in hours
	else 
		ElapsedHumanTime=$( $QueueCommand | grep "$jobID" | awk '{ print $6 }' | tr ':' ' ' | awk '{ print ($3+$2*60.+$1*3600.)/(60.*60.) }' ) #in hours
	fi

	#### COMPUTE THE EFFICIENCY ###
	Efficiency=$( python -c "print float($ElapsedPhysicalTime)/float($ElapsedHumanTime) * 24." )

	echo "Elapsed SCF iterations: ${ElapsedPhysicalTime} iterations."
	echo "Elapsed Human Time: ${ElapsedHumanTime} hours"
	echo "=== EFFICIENCY: $Efficiency iterations/day. ==="
done


#
#squeue | grep 'tjyderri'


