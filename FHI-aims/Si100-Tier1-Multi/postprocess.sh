#!/bin/bash

# set -e

rm energy-vs-kgrid* || true

cd ../Tutorial2/reference/; source tut2.env; cd -

ESiBulk=-7874.14455898 #-7874.14455465 #-7874.14447238 #-7874.14556058 #eV

python_calc () {
    python -c "print ($@)"
}

grep_energy () {
    E=$(grep 'Total energy corrected' <$1 \
	  | gawk '{ printf "%20.8f", $6; }')
    N=$(grep 'Number of atoms' <$1 \
	  | gawk '{print $6; }')
    NSi=$( grep ': Species Si' aims.out | wc -l )
    python_calc "($E)/float($NSi)"
}

grep_cells () {
# just to check
    cells=$( pwd | awk '{ print $1 }' | tr -s "\/" " " | tr -s 'Si-' ' ' | tr -s 'k ' ' ' | tr -s 'x' ' ' | tr -s 'L' ' ' | tr '.' ' ' | awk '{ print $12 }' )
    echo $cells
}

grep_numberatoms () {
    N=$(grep 'Number of atoms' <$1 \
          | gawk '{print $6; }')
    echo $N
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

grep_surface_energy () {
    # calculate surface energy for non-terminated material
    E=$(grep 'Total energy corrected' <$1 \
          | gawk '{ printf "%20.8f", $6; }')
    N=$(grep 'Number of atoms' <$1 \
          | gawk '{print $6; }')
    cells=$( pwd | awk '{ print $1 }' | tr -s "\/" " " | tr -s 'Si-' ' ' | tr -s 'k ' ' ' | tr -s 'x' ' ' | tr -s 'L' ' ' | tr '.' ' ' | awk '{ print $12 }' )
    python_calc "(float($E)-float($N)*($ESiBulk))/(float(2*$cells*$cells))"
}


for dir in `du | awk '{ print $2 }' | grep 'k-'`; do
    pushd $dir > /dev/null
#      echo "# lattice constants in AA, energies in eV" > energy-vs-lattconst.dat
#      for A in [0-9]*.[0-9]*; do
	  # grep --with-filename '*' aims.out || true
	  energy=$( grep_energy aims.out )
          totalenergy=$( grep_total_energy aims.out )
          surfaceenergy=$( grep_surface_energy aims.out )
          numberofcells=$( grep_cells )
          numberofatoms=$( grep_numberatoms aims.out )
	  echo $dir $energy $totalenergy $surfaceenergy $numberofcells $numberofatoms
#      done
#    sort -k1g >> energy-vs-lattconst.dat
    popd > /dev/null
done | sort -k1g >> energy-vs-kgrid.dat 

# Filter the energy list
echo "# Si100. k-point grid= 1 x K x K. Naming convention:<xyz>. Perp to surface=<x>, parallel to surface=<y,z>." > energy-vs-kgrid-filtered.dat
echo "# Blocks SizeXY k-points Vacuum EnergyPerAtom TotalEnergy SurfaceEnergy" >> energy-vs-kgrid-filtered.dat
# Let's build it using two blocks with awk. 
cat energy-vs-kgrid.dat | awk '{ print $1 }' | tr -s "\/" " " | tr -s 'Si-' ' ' | tr -s 'k ' ' ' | tr -s 'x' ' ' | tr -s 'L' ' ' | tr '.' ' ' | awk '{ print $1, $2, $7, $4 }' > Field1.tmp
cat energy-vs-kgrid.dat | awk '{ if($5) {print $2, $3, $4;} else {print "?" }}' > Field2.tmp
paste Field1.tmp Field2.tmp | sort -k1 -k2 -k3 -g >> energy-vs-kgrid-filtered.dat 
rm Field*.tmp
