#!/bin/bash

backup_StH_fixed() {
 DATE=`date +%Y-%m-%d`
 echo "== Saving the StH pulses, not the STH." 
 filename=$DATE-StH-WellControlled-Pulses.tar
 tar --exclude='*.obf' --exclude='core.draco*' -cvf $filename $( list_folders | grep 'StH' )
}

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
 tar -cvf $filename $( find . -name "*.eps" | grep 'SinSquared' )
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
  filename=$DATE-SinSquared-Pulses.tar
  tar --exclude='*.obf' -cvf $filename $( list_folders | grep 'SinSquared' -i )
}

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

remove_td_restart() {
#  find . -mindepth 5 | grep 'lambda' | grep 'obf' >& CanBeRemoved.log
  echo "Warning: THE FOLLOWING FILES WOULD BE REMOVED."
  find . -mindepth 5 -name '*.obf' | grep 'lambda' | tail
  echo "Warning: THE FOLLOWING LINES SHOULD SHOW EMPTY SEARCH!"
  find . -mindepth 5 -name '*.obf' | grep -v 'lambda' | tail
  echo "If you are ready, type find . -mindepth 5 -name '*.obj' #-delete"

}
