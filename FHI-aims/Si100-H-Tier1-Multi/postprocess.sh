#!/bin/bash

# set -e

rm energy-vs-kgrid* || true

cd ../Tutorial2/reference/; source tut2.env; cd -

ESiBulk=-7874.20351239 #eV

python_calc () {
    python -c "print ($@)"
}

grep_energy () {
    E=$(grep 'Total energy corrected' <$1 \
	  | gawk '{ printf "%20.8f", $6; }')
    N=$(grep 'Number of atoms' <$1 \
	  | gawk '{print $6; }')
    NSi=$( grep ': Species Si' <$1 | wc -l )
    python_calc "($E)/float($NSi)"
}

grep_total_energy () {
    E=$(grep 'Total energy corrected' <$1 \
          | gawk '{ printf "%20.8f", $6; }')
#    N=$(grep 'Number of atoms' <$1 \
#          | gawk '{print $6; }')
#    NSi=$( grep ': Species Si' aims.out | wc -l )
#    python_calc "$E"
     echo $E
}

grep_cells () {
# just to check
    cells=$( pwd | awk '{ print $1 }' | tr -s "\/" " " | tr -s 'Si-' ' ' | tr -s 'k ' ' ' | tr -s 'x' ' ' | tr -s 'L' ' ' | tr '.' ' ' | awk '{ print $12 }' )
    echo $cells
}

grep_numberofatoms () {
    grep 'Number of atoms' <$1 \
      | gawk '{print $6; }'
}

for dir in `du | awk '{ print $2 }' | grep 'k-'`; do
    pushd $dir > /dev/null
#      echo "# lattice constants in AA, energies in eV" > energy-vs-lattconst.dat
#      for A in [0-9]*.[0-9]*; do
	  # grep --with-filename '*' aims.out || true
	  energy=$( grep_energy aims.out )
          totalenergy=$( grep_total_energy aims.out )
#          NumberOfAtoms=$( grep_number_of_atoms aims.out )
	  echo ${dir} ${energy} ${totalenergy} 
#      done
#    sort -k1g >> energy-vs-lattconst.dat
    popd > /dev/null
done | sort -k1g >> energy-vs-kgrid.dat 

# Filter the energy list
echo "# Si100. k-point grid= 1 x K x K. Naming convention:<xyz>. Perp to surface=<x>, parallel to surface=<y,z>." > energy-vs-kgrid-filtered.dat
echo "# Layers SizeXY k-points Vacuum EnergyPerAtom TotalEnergy" >> energy-vs-kgrid-filtered.dat
# Let's build it using two blocks with awk. 
cat energy-vs-kgrid.dat | awk '{ print $1 }' | tr -s "\/" " " | tr -s 'Si-' ' ' | tr -s 'k ' ' ' | tr -s 'x' ' ' | tr -s 'L' ' ' | tr '.' ' ' | awk '{ print $1, $2, $7, $4 }' > Field1.tmp
cat energy-vs-kgrid.dat | awk '{ if($2) {print $2, $3;} else {print "?" }}' > Field2.tmp
paste Field1.tmp Field2.tmp | sort -k1 -k2 -k3 -g >> energy-vs-kgrid-filtered.dat 
rm Field*.tmp
