#!/bin/bash

module load git
# module load gnuplot
module load anaconda/2_4.1.1
module list

source $QuantumLaP/02-Slab2D/functions.sh
 
rm *.tmp 

print "Computing: int d^3 r for each available timestep."
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
dt=$( grep 'TimeStep' ../inp | awk '{ print $3 }' )
AUtoSI_time="2.44e-17"
for i in $( du --max-depth=1 | awk '{ print $2 }' | head -n-1 ); 
do
	name=$( cat $i/folders.tmp )
	timestepnumber=$( echo $name | tr "./" " " | awk '{ print $2*1. }'  ) 
        instantAU=$( python_calc " $timestepnumber * $dt " )
	instantSI=$( python_calc " $instantAU * $AUtoSI_time " )
	density=$( cat $i/densities.tmp ); 
	echo "$timestepnumber $instantAU $instantSI $density" >> postCalculationNexcited.tmp 
done

echo "Sorting the final file postCalculationNexcited.log."
echo "# TimeStep \t t (au) \t t (s) \t Ne (au)"
sort -n postCalculationNexcited.tmp >> postCalculationNexcited.log
