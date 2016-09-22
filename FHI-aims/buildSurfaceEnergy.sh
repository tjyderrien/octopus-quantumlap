#!/bin/bash

module load gnuplot
# Purpose: calculate the surface energy with H atoms, and check convergence as function of layers. 

python_calc () {
    python -c "print ($@)"
}

grep_Nideal () {
    # calculates the number of atoms in ideal (non-reconstructed slab)
    AtomPerBlock=python_calc "4*$m*$m"
    blocks=$( awk '{ print $1 }' <$1 )
    python_calc "($AtomPerBlock)*float($blocks)"
}

# Folder names
EnergyFile="energy-vs-kgrid-filtered.dat"
Si100="Si100-Tier1-Multi/$EnergyFile"
Si100H="Si100-H-Tier1-Multi/$EnergyFile"

ESi1=-7874.2035139 #eV
blocks=8 #atom per block

# Just an example: TODO with scripting
ESi=-251971.50823585
ESiH=-252033.97973069
n=1; m=1; 

echo "Nideal..."
Nideal=$( grep_Nideal $Si100 )
echo $Nideal

echo "Surface energy without H"
# Esurf=$( gnuplot -e "print 0.5/(${n}*${m}) * ( ${ESi} - ( ${ESi1} * ${Nideal}  ) )" ) 
#Stekolnikov et al. PRB (2002)

# echo $Esurf

