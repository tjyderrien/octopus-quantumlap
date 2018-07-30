#!/bin/bash
# This script disables the FromScratch of all the octopus runs present below the present folder.
for folder in $( list_folders | grep 'td.general' ); do grep 'FromScratch' $folder/../inp; done
