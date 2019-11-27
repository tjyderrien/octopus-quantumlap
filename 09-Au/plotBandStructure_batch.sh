#!/bin/bash
for folder in $( find | grep 'bandstructure' | grep 'PBE' | sed 's/bandstructure//g' )
do 
	cd $folder; 
	gnuplot '/home/tjyderrien/octopus-slabs/06-Mo/plotBandStructure.gnu'
 	cd -
done
