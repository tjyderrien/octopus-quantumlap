#!/bin/bash
source $QuantumLaPsources/machines.env

CheckComputerName # if EOS: ModuleForEOS 

echo "/!\ Visualization is necessary to obtain the pictures."
# echo "Hence, open visualization server using ./vncRemoteVisualization.sh file."

rm *.mpg *.png

echo "1. Use buildListOfFiles.sh"
echo "2. Correct the file manually."
echo "3. Check plotDensity.py, then run plotDensity.sh."
echo "4. DONT Git push to repository!! Use rsync and GWDG VPN access to get your files. "

echo "Computer name: $ComputerName"
