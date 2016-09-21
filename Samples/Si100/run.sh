#!/bin/bash
#set -e

if [ -e $1 ]; then
	echo "Usage: ./write-geometry.py <Number of Si layers> <Number of cells // surface>"
        echo "Generates FHI-aims input geometry.in file."
        echo "Also generates the input for Octopus."
else 
	./write-geom.py $1 $2 > geometry.in
        ../fhi2octopus.sh geometry.in > geometry.xyz
	echo "You can visualize 'geometry.in' with jmol."
fi
