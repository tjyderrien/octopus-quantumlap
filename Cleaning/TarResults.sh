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
 tar -cvf $filename --exclude='*.obf' $( list_folders | grep 'StH' )
}

backup_eps() {
 DATE=`date +%Y-%m-%d`
 echo "== Saving the EPS files" 
 filename=$DATE-StH-AllEPS.tar
 tar -cvf $filename $( find . -name "*.eps" | grep 'StH' | grep 'eps' )
}
