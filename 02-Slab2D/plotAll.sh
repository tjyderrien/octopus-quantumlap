#!/bin/bash
module load git
module load gnuplot
module load anaconda/2_4.1.1
module list

echo "Building the list of folders can be long..."
for folder in $( du | grep 'td.general' | awk '{ print $2 }' )
do
	echo "Found folder $folder"
	cd $folder
	if [ -e laser ]; then
		gnuplot "${QuantumLaP}/02-Slab2D/plotTimeLaser.gnu"
	fi
	if [ -e energy ]; then
		gnuplot "${QuantumLaP}/02-Slab2D/plotTimeEnergy.gnu"
	fi
	if [ -e total_current ]; then
		gnuplot "${QuantumLaP}/02-Slab2D/plotTimeCurrent.gnu"
	fi
	if [ -e n_ex ]; then
#		gnuplot "${QuantumLaP}/02-Slab2D/plotTimeDensity.gnu"
#		echo "Plotting density..."
                ${QuantumLaP}/02-Slab2D/plotTimeDensityField.py
	fi
	# git add *.eps *.png
	cd - > /dev/null
done

#git commit -m 'Added new time-dependent plots.'
echo "plots have been added to Git. Just commit and push them."

# echo "building a TAR archive also..."
tar -cvf LastEPSfiles.tar $( find | grep 'eps' )
echo "Transfer:"
echo "rsync -avzhe ssh --progress LastEPSfiles.tar thibault@134.76.2.54:/home/thibault/Documents/octopus/02-Slab2D/"

