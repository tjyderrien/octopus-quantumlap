#!/bin/bash

for folder in $( du --max-depth=1 | awk '{ print $2 }' | head -n-1 )
do
	cd $folder
	echo $folder
	$QuantumLaP/02-Slab2D/compareNexcitedWithNeVariations.py
	cd - > /dev/null
done
