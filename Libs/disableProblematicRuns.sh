#!/bin/bash

echo "Usage: ./filterProblematicRuns.sh <filename_with_paths>"
set -e
filename=$1

for folder in $( cat ${filename} )
do
    echo "Treating $folder ..."
    oldpwd=$( pwd )
    cd ${folder}
    recent=$(ls -tr *octopus-td*.log | tail -n1)
    echo "Most recent LOG file: $recent"
    read -s -n 1 -p "Press any key to disable it . . ."
    echo " " > DISABLED && echo "Created file DISABLED in run ${folder}"
#     rm DISABLED && echo "Removed the file located at wrong place."
    cd ${oldpwd}
done

