#!/bin/bash

module load git
module load gnuplot
module load anaconda
module list

source $QuantumLaPsources/Scripts/functions.sh
 
rm *.tmp 

echo "Computing: int d^3 r for each available timestep."
echo ""
echo "Building up list of folders..."
# TotalNumberOfFolders=$( du --max-depth=1 | awk '{ print $2 }' | head -n-1 | wc -l )
TotalNumberOfFolders=$( ls -1 | grep 'td.' | wc -l )
echo "** Total number of folders: $TotalNumberOfFolders"

# if [ ! -e postCalculationNexcited.log ]; then
	counter=0
	for folder in $( ls -1 | grep 'td.' )
	do
		echo $folder
	#	echo "$folder" >> folders.tmp
		previousFolder=$( pwd )
		cd $folder
		#counter=$( python_calc $counter+1 )
		counter=$(( counter + 1 )) #TODO: better perf?
	        echo "** Folder ${counter} / ${TotalNumberOfFolders}."
		echo "$folder" > folders.tmp
		$OctopusSlabsSources/02-Slab2D/compareNexcitedWithNeVariations.py > densities.tmp
		cd $previousFolder > /dev/null
	done

	echo "BREAK HERE, wait for the other"
	while [ $( top -n1 | grep 'tjyderri' | grep 'compareNexcitedWithNeVariations.py' | wc -l ) -ne "0" ]
	do
		sleep 1
		echo "Still $( top -n1 | grep 'tjyderri' | grep 'compareNexcitedWithNeVariations.py' | wc -l ) processes running..."
	done
	echo "All processes finished."

	echo "Collecting the files..."
	dt=$( grep 'TimeStep' ../inp | awk '{ print $3 }' )
	AUtoSI_time="2.418884326505e-17"

	counter=0
	for i in $( ls -1 | grep 'td.' ); 
	do
		counter=$( python_calc $counter+1 )
		echo "** Folder ${counter} / ${TotalNumberOfFolders}."
		name=$( cat $i/folders.tmp )
		timestepnumber=$( echo $name | tr "./" " " | awk '{ print $2*1. }'  ) 
	        instantAU=$( python_calc " $timestepnumber * $dt " )
		instantSI=$( python_calc " $instantAU * $AUtoSI_time " )
		ElectronDensity=$( cat $i/densities.tmp | awk '{ print $1 }' ); 
		HoleDensity=$( cat $i/densities.tmp | awk '{ print $2 }' );
		TotalDensity=$( cat $i/densities.tmp | awk '{ print $3 }'  );
		TotalDensity2=$( cat $i/densities.tmp | awk '{ print $5 }'  );
		echo "$timestepnumber $instantAU $instantSI $ElectronDensity $HoleDensity $TotalDensity $TotalDensity2" >> postCalculationNexcited.tmp 
	done

	echo "Sorting the final file postCalculationNexcited.log."
	echo "# TimeStep \t t (au) \t t (s) \t Ne (au) \t Nh (au) \t Ntot (au) " > postCalculationNexcited.log
	sort -n postCalculationNexcited.tmp >> postCalculationNexcited.log
# fi

echo "Plotting the comparison with excited density..."
gnuplot $OctopusSlabsSources/02-Slab2D/compareNexcitedWithNeVariations.gnu #--persist

