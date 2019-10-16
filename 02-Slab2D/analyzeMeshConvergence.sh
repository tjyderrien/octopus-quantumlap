#!/bin/bash
source $QuantumLaPsources/02-Slab2D/machines.env
CheckComputerName

# ListOfRuns="Si-4x4x6-k-1x1x1-L12-Mesh0.434637008855-Mesh0.188972612546/ Si-4x4x6-k-1x1x1-L12-Mesh0.434637008855-Mesh0.226767135055/ Si-4x4x6-k-1x1x1-L12-Mesh0.434637008855-Mesh0.264561657564/ Si-4x4x6-k-1x1x1-L12-Mesh0.434637008855-Mesh0.302356180073/ Si-4x4x6-k-1x1x1-L12-Mesh0.434637008855-Mesh0.340150702582/ Si-4x4x6-k-1x1x1-L12-Mesh0.434637008855-Mesh0.377945225092/ Si-4x4x6-k-1x1x1-L12-Mesh0.434637008855-Mesh0.415739747601/ Si-4x4x6-k-1x1x1-L12-Mesh0.434637008855-Mesh0.434637008855/"

for runfolder in $0
do
	cd ${runfolder}
	echo "Folder: ${runfolder}"
	grep -H 'rel_ev' ${ResultFile} | tail -n1
	cd - > /dev/null
done > MeshConvergence.dat
echo "Wrote MeshConvergence.dat"
