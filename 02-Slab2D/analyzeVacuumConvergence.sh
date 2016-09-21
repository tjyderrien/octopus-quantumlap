#!/bin/bash
source machines.env
CheckComputerName

ListOfRuns="Si-4x4x6-k-1x1x1-L14-Mesh0.434637008855 Si-4x4x6-k-1x1x1-L12-Mesh0.434637008855 Si-4x4x6-k-1x1x1-L10-Mesh0.434637008855 Si-4x4x6-k-1x1x1-L8-Mesh0.434637008855 Si-4x4x6-k-1x1x1-L6-Mesh0.434637008855"

echo "Creating VacuumConvergence.dat.."
for runfolder in ${ListOfRuns}
do
	cd ${runfolder}
	echo "Folder: ${runfolder}"
	grep -H 'rel_ev' ${ResultFile} | tail -n1
	cd - > /dev/null
done > VacuumConvergence.dat

echo ""

echo "Verification of the potential components convergence."
for runfolder in ${ListOfRuns}
do
	cd ${runfolder}/static/
	pwd
	grep 'Energy' info -A 17
	cd - > /dev/null
done > VacuumConvergenceDetails.dat
ls -lhtr VacuumConverge*.dat

# echo "Check VacuumConvergencDetails.dat"
