#!/bin/bash

#DryRun="--dry-run"
DryRun=""

set -e

if [ -z $2 ]; then
	echo "Usage: ./DuplicateGroundState.sh <Source folder/> <Target folder/>"
	exit
fi

echo "rsync --dry-run -avher $1 $2"
#echo "Add --no-dry-run to confirm the transfer."

if [ $3=="" ]; then
	rsync ${DryRun} -avher $1 $2
else 
	if [ $3=="--no-dry-run" ]; then
		rsync -avher $1 $2
	fi
fi
