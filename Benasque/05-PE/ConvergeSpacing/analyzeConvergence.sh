#!/bin/bash

for folder in $( du --max-depth=1 | awk '{ print $2 }' )
do
	cd $folder
	energy=$( grep 'abs_ev' run.txt | awk '{ print $3 }' | tail -n1 )
	echo "$folder $energy"
	cd - > /dev/null
done

