#!/bin/bash
for folder in $( find . -maxdepth 3 -mindepth 1 | grep 'bandstructure' | sed 's/bandstructure//g' )
do 
	cd $folder; 
	gnuplot "$QuantumLaPsources/Scripts/plotBandStructure.gnu"
 	cd -
done
