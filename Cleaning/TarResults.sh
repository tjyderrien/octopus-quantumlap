#!/bin/bash

backup_STH(){

 echo "Run me from the folder you would like to save." #Save the folder given in argument. "
 DATE=`date +%Y-%m-%d`
 filename=$DATE-STH-Pulses.tar
 
 echo "== Saving the STH pulses, not the StH." 
 echo $DATE
 echo Filename: $filename
 tar -cvf $filename --exclude='*.obf' $( list_folders | grep 'STH' | grep -v 'StH' )

}

backup_StH_fixed() {
 DATE=`date +%Y-%m-%d`
 echo "== Saving the StH pulses, not the STH." 
 filename=$DATE-StH-WellControlled-Pulses.tar
 tar -cvf $filename --exclude='*.obf core.draco*' $( list_folders | grep 'StH' )
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

backup_SinSquared() {
  DATE=`date +%Y-%m-%d`
  echo "== Saving the SinSquared pulses."
  filename=$DATE-SinSquared-Pulses.tar
  tar -cvf $filename --exclude='*.obf' $( list_folders | grep 'SinSquared' -i )
}
