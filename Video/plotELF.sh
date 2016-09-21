#!/bin/bash
echo "/!\ Visualization is necessary to obtain the pictures."
echo "Hence, open visualization server using ./vnc.sh file."
source $QuantumLaP/Video/loadModules.sh
ModuleForEOS
#module load visit
rm ELF*.png ELF*.mpg
visit -cli -s /scratch/tjyderri/octopus/Video/plotELF.py
