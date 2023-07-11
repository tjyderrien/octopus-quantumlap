#!/bin/bash
for folder in $( find . -maxdepth 4 -mindepth 1 | grep 'bandstructure' | sed 's/bandstructure//g' )
do 
	cd $folder; 
	gnuplot "$QuantumLaPsources/octopus_slabs/Scripts/plotBandStructure.gnu"
 	cd -
done
