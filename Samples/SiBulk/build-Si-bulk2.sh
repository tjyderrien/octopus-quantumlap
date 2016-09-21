#!/bin/bash

set -e
# Stop on error
# for A in 5.4; do

A=5.43
echo "Size of Si lattice: $A Angstrums"

# Use this construct for simple calculations . As values
# are replaced verbatim , always put them into "(" , ")".
A4=$( python -c "print ( $A )/4. " )
A2=$( python -c "print ( $A )/2. " )

# Write geometry . in
cat << EOF > geometry.in
# diamond structure with lattice constant $A
lattice_vector 	0.0 $A2 $A2
lattice_vector 	$A2 0.0 $A2
lattice_vector 	$A2 $A2 0.0

atom 		0.0  0.0  0.0 Si
atom 		$A4  $A4  $A4 Si 
EOF
# Write control . in

echo "Visualize the structure with >$ jmol geometry.in &"
