#!/bin/bash

DryRun="--dry-run"
# DryRun=$3

set -e

if [ -z $2 ]; then
	echo "Usage: ./DuplicateGroundState.sh <Source folder/> <Target folder/>"
	exit
fi

echo "rsync ${DryRun} -avhre ssh --progress $1 $2"
#echo "Add --no-dry-run to confirm the transfer."

# if [ $3=="" ]; then
	rsync ${DryRun} -avhre ssh --progress $1 $2
# else 
# 	if [ $3 == "--no-dry-run" ]; then
# 		rsync -avher $1 $2
# 	fi
# fi
