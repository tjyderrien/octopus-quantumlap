#!/bin/bash

##### DATA BACKUPS ######

backup_StH_fixed() {
 DATE=`date +%Y-%m-%d`
 echo "== Saving the StH pulses, not the STH." 
 filename=$DATE-StH-WellControlled-Pulses.tar
 list_folders | grep 'StH' >& backup_StH_list.txt
 tar --exclude='*.obf' --exclude='occs' --exclude="restart" --exclude='core.draco*' -cv -T backup_StH_list.txt -f $filename 
}

backup_video() {
  DATE=`date +%Y-%m-%d`
  echo "== Saving the PNG files present in subfolders..."
  filename=$DATE-Video.tar
  tar -cvf $filename $( find . -name "*.png" | grep 'output_iter' ) 
} 

backup_SinSquared() {
  DATE=`date +%Y-%m-%d`
  echo "== Saving the SinSquared pulses."
  filename=$DATE-SinSquared-Pulses.tar.bz2
  list_folders | grep 'SinSquared' -i >& backup_SinSquared.txt
  tar --exclude='*.obf' --exclude="core" --exclude='occs' --exclude="restart" -cjv -T backup_SinSquared.txt -f $filename
}

backup_SinFourth() {
  DATE=`date +%Y-%m-%d`
  echo "== Saving the SinFourth pulses."
  filename=$DATE-SinFourth-Pulses.tar.bz2
  list_folders | grep 'SinFourth' -i >& backup_SinFourth.txt
  tar --exclude='*.obf' --exclude="core" --exclude='occs' --exclude="restart" -cjv -T backup_SinFourth.txt -f $filename
}

backup_matrixelements() {
  DATE=`date +%Y-%m-%d`
  echo "== Saving the matrix elements."
  filename=$DATE-MatrixElements.tar
  list_folders | grep 'static' >& backup_ME.txt
  ls -1 inp >> backup_ME.txt || echo "Input file was not backuped."
  tar --exclude='*.obf' -cv -T backup_ME.txt -f $filename
}

###### EPS BACKUPS ###########

backup_eps_STH() {
 DATE=`date +%Y-%m-%d`
 echo "== Saving the EPS files" 
 filename=$DATE-STH-AllEPS.tar
 tar -cvf $filename $( find . -name "*.eps" | grep 'STH' | grep -v 'StH' | grep 'eps' )
}

backup_eps_StH() {
 DATE=`date +%Y-%m-%d`
 echo "== Saving the EPS files" 
 filename=$DATE-StH-AllEPS.tar
 tar -cvf $filename $( find . -name "*.eps" | grep 'StH' | grep 'eps' )
}

backup_eps_SinSquared() {
 DATE=`date +%Y-%m-%d`
 echo "== Saving the EPS files for SinSquared runs"
 filename=$DATE-SinSquared-AllEPS.tar 
 tar -cvf $filename $( find . -name "*.eps" | grep 'SinSquare' )
}

backup_eps_SinFourth() {
 DATE=`date +%Y-%m-%d`
 echo "== Saving the EPS files for SinFourth runs"
 filename=$DATE-SinFourth-AllEPS.tar
 tar -cvf $filename $( find . -name "*.eps" | grep 'SinFourth' )
}

#### REDUCE THE NUMBER OF FILES

remove_corefiles() {
  for file in $( find | grep 'core.draco0' )
  do
    rm $file
  done
}

remove_obf_files() {
  # find . -type f | grep '.obf' >& OBFfilesToRemove.log
  find . -name '*.obf' -mindepth 5 -maxdepth 5 #-delete  
}

remove_gs_restart() { 
  echo "Removes all restart information from Ground States."
  echo "Very useful to clean the inodes."
  for folder in $(  find . -type d | grep 'restart/gs' | grep -v 'Si-1x1x1-k-24x24x1-Mesh0.434637008855-Mesh0.434637008855-T0e0' | grep -v 'Si-1x1x1-k-24x24x1-Mesh0.377945225092-Mesh0.377945225092-T0e0' )
  do
    cd $folder
    pwd 
    rm *
    pro-fs
    cd - > /dev/null
  done
}

remove_td_restart() {
#  find . -mindepth 5 | grep 'lambda' | grep 'obf' >& CanBeRemoved.log
  echo "Warning: THE FOLLOWING FILES WOULD BE REMOVED."
  find . -mindepth 5 -name '*.obf' | grep 'lambda' | tail
  echo "Warning: THE FOLLOWING LINES SHOULD SHOW EMPTY SEARCH!"
  find . -mindepth 5 -name '*.obf' | grep -v 'lambda' | tail
  echo "If you are ready, type find . -mindepth 5 -name '*.obf' #-delete"

}
