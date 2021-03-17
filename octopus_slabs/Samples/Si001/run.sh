#!/bin/bash
set -e

if [ -e $1 ]; then
	./write-geom.py $1 > geometry.in
	../fhi2octopus.sh geometry.in > geometry.xyz
else 
	echo "Usage: ./write-geometry.py <Number of Si layers>"
	echo "Generates FHI-aims input geometry.in file."
	echo "Also generates the input for Octopus."
fi
