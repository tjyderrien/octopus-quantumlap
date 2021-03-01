#!/bin/bash
set -e

echo "Loading libraries..."
source $QuantumLaPsources/02-Slab2D/machines.env

echo "Getting computer name..."
CheckComputerName

# ListOfRuns="Si-4x4x6-k-1x1x1-L12-Mesh0.434637008855-Mesh0.188972612546/ Si-4x4x6-k-1x1x1-L12-Mesh0.434637008855-Mesh0.226767135055/ Si-4x4x6-k-1x1x1-L12-Mesh0.434637008855-Mesh0.264561657564/ Si-4x4x6-k-1x1x1-L12-Mesh0.434637008855-Mesh0.302356180073/ Si-4x4x6-k-1x1x1-L12-Mesh0.434637008855-Mesh0.340150702582/ Si-4x4x6-k-1x1x1-L12-Mesh0.434637008855-Mesh0.377945225092/ Si-4x4x6-k-1x1x1-L12-Mesh0.434637008855-Mesh0.415739747601/ Si-4x4x6-k-1x1x1-L12-Mesh0.434637008855-Mesh0.434637008855/"

echo "Building up list of folders..."
ListOfRuns=$( find . -type d -maxdepth 2 -mindepth 2 )

echo $ListOfRuns

for runfolder in $ListOfRuns
do
	cd ${runfolder}
	echo "Folder: ${runfolder}"
        ResultFile=$( ls -tr -1 *.log | tail -n1 )
        if [ -a $ResultFile  ]; then
          echo "Treating file $ResultFile"
	  grep -H 'etot' ${ResultFile} | tail -n1 | awk '{ print $4 }' 
          grep -H 'Fermi' -i ${ResultFile}
        else 
          echo "Didnt find anything in folder ${runfolder}"
        fi
	cd - > /dev/null
done > MeshConvergence.dat
echo "Wrote MeshConvergence.dat"
