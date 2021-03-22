#!/bin/bash

# Executed by cron tabs
# crontab -e @weekly ./RegularBackup.sh

# set -e

source $QuantumLaPsources/octopus_slabs/Cleaning/TarResults.sh
source $QuantumLaPsources/02-Slab2D/functions.sh
cd $QuantumLaPruns

for folder in $( find . -mindepth 1 -maxdepth 1 -type d | grep -v 'backuped')
do
	cd $folder
        pwd
        backup_matrixelements
	backup_SinFourth
	backup_SinSquared
	backup_StH_fixed
	backup_video
	cd - > /dev/null
done
