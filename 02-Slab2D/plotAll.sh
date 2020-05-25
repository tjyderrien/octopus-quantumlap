#!/bin/bash
module purge
# module load git
# module load gnuplot #WARNING: GNUplot is now loaded through machines.env.
#module load anaconda/2/5.1 #we shall use machines.env library load for this.
source $QuantumLaPsources/machines.env
CheckComputerName #this loads the machine dependent libraries from machines.env.

echo "** Info: UseScreen? [${UseScreen}]"

source $QuantumLaPsources/02-Slab2D/functions.sh
# module load plgrid/tools/gnuplot/5.2.4 #temporary fix
module list
# set -e

echo "Current folder: $( pwd )"

export OctopusSlabsSources=${QuantumLaPsources}
echo "** Debug: OctopusSlabsSources=${OctopusSlabsSources}"

# export ListOfRuns="lambda1600e-9-tau30E-15-field1E8-CEP0.0-StH2 lambda2200e-9-tau30E-15-field1E8-CEP0.0-StH2 lambda3200e-9-tau30E-15-field1E8-CEP0.0-StH2 lambda1600e-9-tau30E-15-field1E9-CEP0.0-StH2 lambda1600e-9-tau30E-15-field5E9-CEP0.0-StH2 lambda1600e-9-tau30E-15-field1E10-CEP0.0-StH2 lambda2200e-9-tau30E-15-field1E9-CEP0.0-StH2 lambda2200e-9-tau30E-15-field5E9-CEP0.0-StH2 lambda2200e-9-tau30E-15-field1E10-CEP0.0-StH2 lambda3200e-9-tau30E-15-field1E9-CEP0.0-StH2 lambda3200e-9-tau30E-15-field5E9-CEP0.0-StH2 lambda3200e-9-tau30E-15-field1E10-CEP0.0-StH2"

export ListOfRuns="lambda800e-9-tau20E-15-field9E8-CEP0.0-StH2 lambda800e-9-tau20E-15-field9.1E8-CEP0.0-StH2 lambda800e-9-tau20E-15-field9.2E8-CEP0.0-StH2 lambda800e-9-tau20E-15-field9.3E8-CEP0.0-StH2 lambda800e-9-tau20E-15-field9.4E8-CEP0.0-StH2 lambda800e-9-tau20E-15-field9.6E8-CEP0.0-StH2 lambda800e-9-tau20E-15-field9.7E8-CEP0.0-StH2 lambda800e-9-tau20E-15-field9.8E8-CEP0.0-StH2 lambda800e-9-tau20E-15-field9.9E8-CEP0.0-StH2 lambda800e-9-tau20E-15-field1.1E9-CEP0.0-StH2 lambda800e-9-tau20E-15-field1.2E9-CEP0.0-StH2 lambda800e-9-tau20E-15-field1.3E9-CEP0.0-StH2 lambda800e-9-tau20E-15-field1.4E9-CEP0.0-StH2 lambda800e-9-tau20E-15-field1.6E9-CEP0.0-StH2"


plotBasicOutput () {

echo "Building the list of folders can be long..."
for folder in $( list_folders | grep 'td.general' )
do
	echo "Found folder $folder"
	cd $folder
	if [ -e laser ]; then
#		gnuplot "${OctopusSlabsSources}/02-Slab2D/plotTimeDensity.gnu"
		echo "** Fixing bad output octopus problems..."
		${QuantumLaPsources}/fixOutputFile.sh laser
        fi
	if [ -e laser ]; then
		echo "** Plotting laser..."
                if [[ $UseScreen ]]; then
                echo " *** Warning: screen is employed."
		screen -d -m gnuplot "${QuantumLaPsources}/02-Slab2D/plotTimeLaser.gnu"
                screen -d -m gnuplot "${QuantumLaPsources}/03-Bulk3D/plotPolarization.gnu"
                else
                gnuplot "${QuantumLaPsources}/02-Slab2D/plotTimeLaser.gnu"
                gnuplot "${QuantumLaPsources}/03-Bulk3D/plotPolarization.gnu"
                fi
	fi
	if [ -e energy ]; then
		echo "** Plotting energy..."
		if [[ $UseScreen ]]; then
                echo " *** Warning: screen is employed."
                screen -d -m gnuplot "${QuantumLaPsources}/02-Slab2D/plotTimeEnergy.gnu"
                else
                gnuplot "${QuantumLaPsources}/02-Slab2D/plotTimeEnergy.gnu"
                fi
	fi
	if [ -e total_current ]; then
		echo "** Plotting current..."
		if [[ $UseScreen ]]; then
                echo " *** Warning: screen is employed."
                screen -d -m gnuplot "${QuantumLaPsources}/02-Slab2D/plotTimeCurrent.gnu"
                else
                gnuplot "${QuantumLaPsources}/02-Slab2D/plotTimeCurrent.gnu"
                fi
	fi
        if [ -e n_ex ]; then
		echo "Plotting density..."
                if [[ $UseScreen ]]; then
                echo " *** Warning: screen is employed."
                screen -d -m gnuplot "${QuantumLaPsources}/02-Slab2D/plotTimeDensity.gnu"
                else
                gnuplot "${QuantumLaPsources}/02-Slab2D/plotTimeDensity.gnu"
                fi
	fi

        if [ -e dielectric_function ]; then
                echo "Plotting dielectric permittivity..."
                if [[ $UseScreen ]]; then
                echo "*** Warning: screen is employed."
                screen -d -m gnuplot "${QuantumLaPsources}/02-Slab2D/plotDielectricPermittivity.gnu"
                else
                gnuplot "${QuantumLaPsources}/02-Slab2D/plotDielectricPermittivity.gnu"
                fi
        fi
	echo ""
	echo "** NEXT FOLDER"
	echo ""
	# git add *.eps *.png
	cd - >& /dev/null
done

# for each results folder where the code was started. #Assume that ./td.general was created. 
}

plotExcitation () {
for folder in $( list_folders | grep 'td.general' )
do
	# HERE, the whole thing should be executed in background for each folder (1 thread / folder).
 	# then, processes may work on a different run, and hence produce all results at the same time
 	# without interfering.
	echo "Found folder $folder"
	prevFolder=$( pwd )
        cd $folder
	# screen -d -m ${OctopusSlabsSources}/02-Slab2D/compareNexcitedWithNeVariations.sh

        if [ -e n_ex ]; then
          echo "** Plotting the TD-ExcitedDensity results... $(pwd)" 
          if [[ $UseScreen ]]; then
            screen -d -m ${QuantumLaPsources}/02-Slab2D/plotTimeKeldysh.py auto 0 0 >& plotTimeKeldysh.log
          else
            ${QuantumLaPsources}/02-Slab2D/plotTimeKeldysh.py auto 0 0 # >& plotTimeKeldysh.log
          fi
#        screen -d -m ${OctopusSlabsSources}/02-Slab2D/plotTimeKeldysh.py auto 1 0
	  echo ""
        fi

 	# echo "** Generating list of TD-visual outputs..."
 	# TODO : $OctopusSlabsSources/Video/buildFileList.sh #generated listoffiles.txt
        # echo "td.0000000/" > listoffiles.txt
 	# echo "td.0007770/" >> listoffiles.txt
        # screen -d -m $OctopusSlabsSources/Video/plotELF.sh

	cd $prevFolder > /dev/null
done
}

plotHHG () {
# plot results for HHG for SquaredSinus pulses
for folder in $( list_folders | grep 'td.general' )
do
        # HERE, the whole thing should be executed in background for each folder (1 thread / folder).
        # then, processes may work on a different run, and hence produce all results at the same time
        # without interfering.
        echo "Found folder $folder"
        prevFolder=$( pwd )
        cd $folder
        echo "** Building files for post-processed density..."
        # screen -d -m ${OctopusSlabsSources}/02-Slab2D/compareNexcitedWithNeVariations.sh
        if [ -e total_current ]; then
          echo "** Plotting the TD-HHG results... $(pwd)" 
          if [[ $UseScreen ]]; then
            screen -d -m ${QuantumLaPsources}/02-Slab2D/plotTimeHHG.py auto 0 0
          else
            ${QuantumLaPsources}/02-Slab2D/plotTimeHHG.py auto 0 0
          fi
        fi

#        screen -d -m ${OctopusSlabsSources}/02-Slab2D/plotTimeHHG.py auto 1 0
        echo ""

        # echo "** Generating list of TD-visual outputs..."
        # TODO : $OctopusSlabsSources/Video/buildFileList.sh #generated listoffiles.txt
        # echo "td.0000000/" > listoffiles.txt
        # echo "td.0007770/" >> listoffiles.txt
        # screen -d -m $OctopusSlabsSources/Video/plotELF.sh

        cd $prevFolder > /dev/null
done

}

plotLinearAbsorption () {
for folder in $( list_folders | grep 'td.general' | grep 'field0E0' | grep 'Probe' ) #check only the places where we do have a pump-probe
do
        echo "Found folder $folder"
        prevFolder=$( pwd )
        cd $folder
        if [ -e total_current ]; then
          if [ -e gauge_field ]; then
            echo "** Plotting the Linear Response results... $(pwd)" 
            if [[ $UseScreen ]]; then
              screen -d -m ${QuantumLaPsources}/02-Slab2D/plotLinearResponse.py
            else
              ${QuantumLaPsources}/02-Slab2D/plotLinearResponse.py
            fi
          fi
        fi

        echo ""

        cd $prevFolder > /dev/null
done


}
# Cleans the octopus files and prepare simple plots for all. 
# plotBasicOutput

# Prepare the results for the list specified at the beginning of the file
# Basically mid-IR experiments
# for run in $ListOfRuns
# do
# LocalPath=$( pwd )
# cd $run


echo "Plotting basic results..."
plotBasicOutput #cleans up laser files also. 
echo "Plotting Excitation results..."
plotExcitation  #Prepare all plots for n_ex analysis.
echo "Plotting harmonic generation results..."
# plotHHG         #Generate the CSV files with HHG spectra. 
echo "** Plotting linear absorption for folders where it makes sense. $(pwd)..."
plotLinearAbsorption
echo " Have fun with Database.csv and plotHHGspectra.py after this."

# for folder in $( list_folders | grep -v 'lambda' ) 
# do
#  prevFolder=$( pwd )
#  cd $folder

# echo "** Plotting nice figure for post-processed density for $(pwd)..."
${QuantumLaPsources}/02-Slab2D/plotFinalQuantities.py

#   cd $prevFolder > /dev/null
# done 

# git commit -m 'Added new time-dependent plots.'
# echo "plots have been added to Git. Just commit and push them."


echo "Once you are done: build a TAR archive with results."
echo 'tar -cvf LastEPSfiles.tar $( find | grep eps )'
echo "For transfer, use:"
echo "rsync -avzhe ssh --progress LastEPSfiles.tar thibault@134.76.2.54:/home/thibault/Documents/octopus/02-Slab2D/"

