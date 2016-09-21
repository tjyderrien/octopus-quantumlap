#!/bin/bash

set -e

# meshList="0.7 1 1.1 1.2 1.3 1.4 1.5 2 3"
radiusList="10 11 12 13 14 15 20"

for radius in $radiusList
do
	dirName="Sphere$radius"
	echo "Creating directory $dirName"
	mkdir -p $dirName
	echo "Going there."
	cd $dirName
	# For each mesh size
	echo "Copying the submission and the source files."
	cp ../../inp-gs ./inp
	cp ../../submit_mpi.sge .
	# Replace in the inp file
	sed -i.bak s/SphereRadius/$radius/g inp
	qsub submit_mpi.sge
	# 
	cd -
done
