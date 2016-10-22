#!/bin/bash

module load git
module load gnuplot
module load anaconda/2_4.1.1
module list
 
rm *.tmp 

for folder in $( du --max-depth=1 | awk '{ print $2 }' | head -n-1 )
do
	echo $folder
#	echo "$folder" >> folders.tmp
	cd $folder
	echo "$folder" > folders.tmp
	$QuantumLaP/02-Slab2D/compareNexcitedWithNeVariations.py > densities.tmp
	cd - > /dev/null
done

echo "Collecting the files..."

for i in $( du --max-depth=1 | awk '{ print $2 }' | head -n-1 ); 
do
	name=$( cat $i/folders.tmp )
	density=$( cat $i/densities.tmp ); 
	echo "$name $density" >> postCalculationNexcited.tmp 
done

echo "Sorting the final file postCalculationNexcited.log."
sort -n postCalculationNexcited.tmp > postCalculationNexcited.log
