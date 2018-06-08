#!/bin/bash

echo "Run me from the folder you would like to save." #Save the folder given in argument. "
DATE=`date +%Y-%m-%d`

tar -cvf $DATE-STH-Pulses.tar $( find | grep 'STH' | grep -v 'StH' | grep 'restart' -v )
