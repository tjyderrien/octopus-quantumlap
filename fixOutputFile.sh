#!/bin/bash

# Cleans the laser file from octopus in case fields are overlapping.
# if [ -e $i.tmp ]; then
#	echo "$1.tmp already exists." >&2
#	exit 1
# fi

if [ -z $1 ]; then
	echo "Usage: script.sh <File to be fixed>." >&2
	echo "Description: clean a file with overlapping fields of floats."
	exit 1
fi


# mv $1 $1.tmp
perl -pe 's/(?<!e)(-)/ $1/g' $1 > $1.tmp
