#!/bin/bash

echo "Usage: ./filterProblematicRuns.sh <filename_with_paths>"
set -e
filename=$1

for folder in $( cat ${filename} )
do
    oldpwd=$( pwd )
    cd ${folder}; cd ..
    recent=$(ls -tr *octopus-td*.log | tail -n1)
    echo "Most recent LOG file: $recent"
    echo " " > DISABLED && echo "Created file DISABLED in run ${folder}"
    cd ${oldpwd}
done

