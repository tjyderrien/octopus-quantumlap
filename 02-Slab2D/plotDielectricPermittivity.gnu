#!gnuplot

reset

set output 'DielectricPermittivityX.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'

set xlabel 'Frequency {/Symbol w} (eV)'
set ylabel 'Permittivity {/Symbol e}'

Hartree_to_eV=27.21138602

file = "dielectric_function"

p file u ($1*Hartree_to_eV):2 w l lt 1 lc 1 t 'Re({/Symbol e}_x)', \
  file u ($1*Hartree_to_eV):3 w l lt 2 lc 2 t 'Im({/Symbol e}_x)' 

#
#  file u ($1*Hartree_to_eV):4 w l lt 1 lc 2 t 'Re({/Symbol e}_y)', \
#  file u ($1*Hartree_to_eV):5 w l lt 2 lc 2 t 'Im({/Symbol e}_y)', \
#  file u ($1*Hartree_to_eV):6 w l lt 1 lc 3 t 'Re({/Symbol e}_z)', \
#  file u ($1*Hartree_to_eV):7 w l lt 2 lc 3 t 'Im({/Symbol e}_z)'
