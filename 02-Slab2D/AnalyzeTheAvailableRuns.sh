#!/bin/bash
# echo "Analysis of available simulations."
# echo "Usage: ./Analyze.sh <Number oflines per output>, <Depth To analyze>"

DepthToAnalyze=$2
NumberLines=$1

DebugMode=0

source ${OctopusSlabsSources}/02-Slab2D/functions.sh
source ${OctopusSlabsSources}/02-Slab2D/machines.env
CheckComputerName

if [ -z ${DepthToAnalyze} ]; then
	echo "Analysis of available simulations."
	echo "Usage: ./Analyze.sh <Number oflines per output>, <Depth To analyze>"
	exit
fi

FindAvailableJobs() {
	find | grep "$ResultFile"
}

FilterTheMatchingFolder() {
	for i in $( list_folders )
	do
		cd $i > /dev/null 
		if [ ${DebugMode} == 1 ]; then
			echo "Path: $i"
		fi
		# Machine dependent job names... TODO: correct in submit.sbatch and submit.pbs
		case "${ComputerName}" in 
			${Drako}*)
				MostRecentFile=$( ls -lhtr octopus*.log 2> /dev/null | tail -n1 | awk '{ print $9 }' )
			;;

			${EOS}*)
				MostRecentFile=$( ls -ltrh run.txt 2> /dev/null | tail -n1 | awk '{ print $9 }'  )
			;; 

			*)
				echo "Computer name error. Name = ${ComputerName}. This computer is not taken into account."
			;;
		esac
		
		if [ ${DebugMode} == 1 ]; then
			echo "Most recent file = ${MostRecentFile}."
		fi

		if [ ${MostRecentFile} ]; then 
			echo "$i/$MostRecentFile"
			ls "${MostRecentFile}" -lhtr
			# grep '**' $MostRecentFile
			tail -n${NumberLines} $MostRecentFile
			echo "*********"
#		else 
#			echo "MostRecentFile was undefined."
		fi
		cd - > /dev/null
	done
}

# for i in $( FindAvailableJobs )
# do
#	echo $i
	FilterTheMatchingFolder	
#	echo ""  # >> SuccessList.dat
#	grep $i | grep ''
# done

# echo "Number of occurences."
FindAvailableJobs | wc -l
# echo "See SuccessList.dat and ErrorList.dat."
