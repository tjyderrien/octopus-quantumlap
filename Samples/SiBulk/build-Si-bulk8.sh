#!/bin/bash

set -e
# Stop on error
# for A in 5.4; do

A=5.43
echo "Size of Si lattice: $A Angstrums"

# Use this construct for simple calculations . As values
# are replaced verbatim , always put them into "(" , ")".
A14=$( python -c "print ( $A )/4. " )
A24=$( python -c "print ( $A )/2. " )
A34=$( python -c "print ( 3. * $A ) / 4. ")

# Write geometry . in
cat << EOF > geometry.in
atom 0.0  0.0  0.0 Si
atom $A14 $A14 $A14 Si
atom $A34 $A34 $A14 Si
atom $A14 $A34 $A34 Si
atom $A34 $A14 $A34 Si
atom 0.0  $A24 $A24 Si
atom $A24 0.0  $A24 Si
atom $A24  $A24 0.0  Si 

lattice_vector $A 0.0 0.0
lattice_vector 0.0 $A 0.0
lattice_vector 0.0 0.0 $A

EOF
# Write control . in

# echo "Visualize the structure with >$ jmol geometry.in &"
