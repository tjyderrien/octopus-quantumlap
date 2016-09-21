#!/bin/bash
source machines.env
CheckComputerName
if [ ${ComputerName} != ${EOS} ]; then 
	echo "This script is meant to work on EOS."
	exit
fi
for i in $( ./WhoAreTheRunningJobs.sh ); 
do 
	echo $i; 
	echo "| ITERATION STEP: "
	tail -n20 $i | grep 'Elapsed time for SCF step' | more;
	echo "| Energy of latest step"
	tail -n100 $i | grep 'rel_ev'
	echo "" 
done
