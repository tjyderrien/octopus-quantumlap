#!/bin/bash

source machines.env
CheckComputerName

echo "Current machine: ${ComputerName}"

echo "RUNNING JOBS"

GetJobID() {
        ${QueueCommand} | grep 'tjyderri' | grep "${RunningMark}" | awk '{ print $1 }'
}

# ${QueueCommand}
for i in $( GetJobID )
do
	echo $i
	for j in $( find | grep "${i}") #| grep "${ResultFile}"
	do
		dirname $j
	done
done
