#!/bin/bash

# Turns output of FHI geometry to Octopus
#cat $1 |  grep -v 'lattice' | wc -l ;

# octopus strictly works with atomic units now. Conversion will be done here, then. 

if [ -z $1 ]; then
	echo "Usage: ./script.sh <geometry.in>"
	exit
fi

# WARNING: Although octopus announced that working in Angstroms in over in next versions, it is actually WRONG stateent for the XYZ files. If using XYZ input, then IT MUST BE IN ANSTROMS. No need to convert then. 
# cat $1 | grep -v 'lattice' | grep -v -e '^$' | awk '{ Bohr=0.52917721067 } { print G$5 " " $2/Bohr " " $3/Bohr " " $4/Bohr  }'
cat $1 | grep -v 'lattice' | grep -v -e '^$' | awk '{ print G$5 " " $2 " " $3 " " $4  }'
