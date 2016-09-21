#!/bin/bash

set -e

rm energy-vs-kgrid.dat

cd ../Tutorial2/reference/; source tut2.env; cd -

python_calc () {
    python -c "print ($@)"
}

grep_energy () {
    E=$(grep 'Total energy corrected' <$1 \
	  | gawk '{printf "%20.8f", $6; }')
    N=$(grep 'Number of atoms' <$1 \
	  | gawk '{print $6; }')
    python_calc "($E)/float($N)"
}

for dir in `du | awk '{ print $2 }' | grep 'k-'`; do
    pushd $dir > /dev/null
      echo "# lattice constants in AA, energies in eV" > energy-vs-lattconst.dat
#      for A in [0-9]*.[0-9]*; do
	  # grep --with-filename '*' aims.out || true
	  energy=$(grep_energy aims.out)
	  echo $dir $energy
#      done
#    sort -k1g >> energy-vs-lattconst.dat
    popd > /dev/null
done | sort -k1g >> energy-vs-kgrid.dat

# Filter the energy list
echo "# Layers k-points Energy" > energy-vs-kgrid-filtered.dat
cat energy-vs-kgrid.dat | tr "\/" " " | tr 'Si-'  ' ' | tr 'k ' ' ' | tr 'x' ' ' | awk '{ print $4, $7, $8 }' >> energy-vs-kgrid-filtered.dat
