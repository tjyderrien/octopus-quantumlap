#!/bin/bash
echo "Baricenter on Z (depth) axis is in Angstr."
./write-geom.py $1 $2 | grep 'lattice_vector' -v | awk 'BEGIN {Sum=0.} { Sum=Sum+$4 } END {print Sum; } '
echo "Baricenter on X axis is in Angstr."
./write-geom.py $1 $2 | grep 'lattice_vector' -v | awk 'BEGIN {Sum=0.} { Sum=Sum+$2 } END {print Sum; } '
echo "Baricenter on Y axis is in Angstr."
./write-geom.py $1 $2 | grep 'lattice_vector' -v | awk 'BEGIN {Sum=0.} { Sum=Sum+$3 } END {print Sum; } '
