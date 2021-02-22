#!/bin/bash
source $QuantumLaPsources/machines.env
LoadEnvLibraries
source $QuantumLaPsources/Video/loadModules.sh
# $RunVisit 
visit -cli -s $QuantumLaPsources/Video/plotDensityDifference3d.py
