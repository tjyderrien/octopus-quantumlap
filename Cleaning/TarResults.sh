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
 echo "== Saving the StH pulses, not the STH." 
 filename=$DATE-StH-Pulses.tar
 tar -cvf $filename $( list_folders | grep 'StH' )
}
