#!/bin/bash

echo "List of available runs."
echo "Convergence state."
for i in $( find | grep 'Si-4x4x6-k-1x1x1-L' | grep 'mesh' -i | grep 'run.txt'  )
do
	echo $i
	grep 'rel_ev' $i | tail -n1
done
