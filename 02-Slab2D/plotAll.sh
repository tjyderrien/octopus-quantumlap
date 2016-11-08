#!/bin/bash
module load git
module load gnuplot
module load anaconda
module list

source $QuantumLaP/02-Slab2D/functions.sh
# set -e

echo "Current folder: $( pwd )"

echo "Building the list of folders can be long..."
for folder in $( list_folders | grep 'td.general' )
do
	echo "Found folder $folder"
	cd $folder
	if [ -e laser ]; then
		echo "** Plotting laser..."
		gnuplot "${QuantumLaP}/02-Slab2D/plotTimeLaser.gnu"
	fi
	if [ -e energy ]; then
		echo "** Plotting energy..."
		gnuplot "${QuantumLaP}/02-Slab2D/plotTimeEnergy.gnu"
	fi
	if [ -e total_current ]; then
		echo "** Plotting current..."
		gnuplot "${QuantumLaP}/02-Slab2D/plotTimeCurrent.gnu"
	fi
	if [ -e n_ex ]; then
#		gnuplot "${QuantumLaP}/02-Slab2D/plotTimeDensity.gnu"
		echo "Fixing bad output octopus problems..."
		${QuantumLaP}/fixOutputFile.sh laser
		echo "Plotting density..."
                ${QuantumLaP}/02-Slab2D/plotTimeDensityField.py
	fi
	echo ""
	echo "** NEXT FOLDER PLEASE :)"
	echo ""
	# git add *.eps *.png
	cd - > /dev/null
done

for folder in $( list_folders | grep 'output_iter' | grep -v 'td.' )
do
	echo "Found folder $folder"
	prevFolder=$( pwd )
        cd $folder
	echo "** Building files for post-processed density..."
#	${QuantumLaP}/02-Slab2D/compareNexcitedWithNeVariations.sh
	echo "** Plotting nice figure for post-processed density for $(pwd)..."
	${QuantumLaP}/02-Slab2D/plotTimeKeldysh.py
	echo ""
	
	cd $prevFolder > /dev/null
done

# git commit -m 'Added new time-dependent plots.'
# echo "plots have been added to Git. Just commit and push them."
echo "Don't use Git for transfering binary data!"

echo "building a TAR archive with results..."
echo 'tar -cvf LastEPSfiles.tar $( find | grep eps )'
echo "Transfer:"
echo "rsync -avzhe ssh --progress LastEPSfiles.tar thibault@134.76.2.54:/home/thibault/Documents/octopus/02-Slab2D/"

