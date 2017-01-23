#!/bin/bash
echo "/!\ Visualization is necessary to obtain the pictures."
echo "Hence, open visualization server using ./vnc.sh file."
source $QuantumLaP/Video/loadModules.sh
CheckComputerName # if EOS: ModuleForEOS 

#module load visit
rm *.mpg *.png

echo "1. Use buildListOfFiles.sh"
echo "2. Correct the file manually."
echo "3. Check plotDensity.py, then run plotDensity.sh."
echo "4. DONT Git push to repository!! Use rsync and GWDG VPN access to get your files. "

echo "Computer name: $ComputerName"


$RunVisit -cli -s $QuantumLaP/Video/plotDensityDifference3d.py
