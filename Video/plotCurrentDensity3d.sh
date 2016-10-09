#!/bin/bash
echo "/!\ Visualization is necessary to obtain the pictures."
echo "Hence, open visualization server using ./vnc.sh file."
source $QuantumLaP/Video/loadModules.sh
ModuleForEOS
=======
#module load visit
rm *.mpg *.png

echo "1. Use buildListOfFiles.sh"
echo "2. Correct the file manually."
echo "3. Check plotDensity.py, then run plotDensity.sh."
echo "4. Git push to repository :) "
visit -cli -s /scratch/tjyderri/octopus/Video/plotCurrentDensity3d.py
