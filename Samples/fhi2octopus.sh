#!/bin/bash

# Turns output of FHI geometry to Octopus
#cat $1 |  grep -v 'lattice' | wc -l ;

# octopus strictly works with atomic units now. Conversion will be done here, then. 

if [ -z $1 ]; then
	echo "Usage: ./script.sh <geometry.in>"
	exit
fi

cat $1 | grep -v 'lattice' | grep -v -e '^$' | awk '{ Bohr=0.52917721067 } { print $5 " " $2/Bohr " " $3/Bohr " " $4/Bohr  }'

