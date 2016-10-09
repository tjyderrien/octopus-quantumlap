#!/bin/bash

# Turns output of FHI geometry to Octopus
#cat $1 |  grep -v 'lattice' | wc -l ;

if [ -z $1 ]; then
	echo "Usage: ./script.sh <geometry.in>"
	exit
fi
cat $1 |  grep -v 'lattice' | awk '{ print $5 " " $2 " " $3 " " $4  }'

