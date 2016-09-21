#!/bin/bash

# Turns output of FHI geometry to Octopus
#cat $1 |  grep -v 'lattice' | wc -l ;
cat $1 |  grep -v 'lattice' | awk '{ print $5 " " $2 " " $3 " " $4  }'

