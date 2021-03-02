#!/bin/bash

# T. J.-Y. Derrien, HiLASE Centre, Institute of Physics <derrien@fzu.cz>
echo "** Generate list of files to be transfered. From biggest to smallest."
for files in $( ls -S1 *AllEPS.tar* )
do
	echo $files
done >& ListOfBackupedFiles.log

echo "** Transfering each file *AllEPS.tar* to MetaCentrum."
for files in $( cat ListOfBackupedFiles.log )
do
    echo "Transfering: $files"
    echo "Password for MetaCentrum?"
    rsync -avzhre ssh --progress $files derrit01@ssh.du4.cesnet.cz:/home/derrit01/VO_metacentrum-tape_tape-archive/migrated_du2/Si/Laptop-AllEPS-Archives/
    rm $files -i
done
