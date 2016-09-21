#!/bin/bash

set -e

meshList="0.7 1 1.1 1.2 1.3 1.4 1.5 2 3"

for mesh in $meshList
do
	dirName="Mesh$mesh"
	echo "Creating directory $dirName"
	mkdir -p $dirName
	echo "Going there."
	cd $dirName
	# For each mesh size
	echo "Copying the submission and the source files."
	cp ../../inp-gs ./inp
	cp ../../submit_mpi.sge .
	# Replace in the inp file
	sed -i.bak s/MeshSize/$mesh/g inp
	qsub submit_mpi.sge
	# 
	cd -
done
