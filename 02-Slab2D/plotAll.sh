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
#		echo "Plotting density..."
#               ${QuantumLaP}/02-Slab2D/plotTimeKeldysh.py
	fi
	echo ""
	echo "** NEXT FOLDER PLEASE :)"
	echo ""
	# git add *.eps *.png
	cd - > /dev/null
done

for folder in $( list_folders | grep 'output_iter' | grep -v 'td.' )
do
	# HERE, the whole thing should be executed in background for each folder (1 thread / folder).
 	# then, processes may work on a different run, and hence produce all results at the same time
 	# without interfering.
	echo "Found folder $folder"
	prevFolder=$( pwd )
        cd $folder
	echo "** Building files for post-processed density..."
        echo "** ENABLE ONLY IF NECESSARY. Very long.. (1h)"
	# screen -d -m ${QuantumLaP}/02-Slab2D/compareNexcitedWithNeVariations.sh
	echo "** Plotting nice figure for post-processed density for $(pwd)..."
	screen -d -m ${QuantumLaP}/02-Slab2D/plotTimeKeldysh.py
	echo ""

 	# echo "** Generating list of TD-visual outputs..."
 	# TODO : $QuantumLaP/Video/buildFileList.sh #generated listoffiles.txt
        # echo "td.0000000/" > listoffiles.txt
 	# echo "td.0007770/" >> listoffiles.txt
        # screen -d -m $QuantumLaP/Video/plotELF.sh

	cd $prevFolder > /dev/null
done

# for folder in $( list_folders | grep -v 'lambda' ) 
# do
#  prevFolder=$( pwd )
#  cd $folder
  echo "** Plotting nice figure for post-processed density for $(pwd)..."
  screen -d -m ${QuantumLaP}/02-Slab2D/plotFinalQuantities.py 0 0
  screen -d -m ${QuantumLaP}/02-Slab2D/plotFinalQuantities.py 0 1
  screen -d -m ${QuantumLaP}/02-Slab2D/plotFinalQuantities.py 1 0
  screen -d -m ${QuantumLaP}/02-Slab2D/plotFinalQuantities.py 1 1
#   cd $prevFolder > /dev/null
# done 

# git commit -m 'Added new time-dependent plots.'
# echo "plots have been added to Git. Just commit and push them."
echo "Don't use Git for transfering binary data!"

echo "building a TAR archive with results..."
echo 'tar -cvf LastEPSfiles.tar $( find | grep eps )'
echo "Transfer:"
echo "rsync -avzhe ssh --progress LastEPSfiles.tar thibault@134.76.2.54:/home/thibault/Documents/octopus/02-Slab2D/"

