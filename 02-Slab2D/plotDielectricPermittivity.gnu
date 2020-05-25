#!gnuplot

reset

set output 'DielectricPermittivityX.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'

set xlabel 'Frequency {/Symbol w} (eV)'
set ylabel 'Permittivity {/Symbol e}'

Hartree_to_eV=27.21138602
epsilon_0 = 8.854187817620389e-12
hbar=1.0545718001391127e-34
ec  =1.6021766208e-19
# set xrange [:]

file = "dielectric_function"

# dielectric=1+conductivity(real,imag)*{0,1}/omega/epsilon_0
conductivity(ReEps,ImEps,omega)={0,1}*epsilon_0*omega*({1,0}-{1,0}*ReEps-{0,1}*ImEps) #SI
# conductivity(ReEps,ImEps,omega)={0,1}*omega/(4.*pi)*({1,0}-{1,0}*ReEps-{0,1}*ImEps) #AU

p file u ($1*Hartree_to_eV):2 w l lt 1 lc "red" t 'Re({/Symbol e}_x)', \
  file u ($1*Hartree_to_eV):3 w l lt 2 lc "blue" t 'Im({/Symbol e}_x)' 

set output 'DielectricPermittivityX-Zoom.eps'
set xrange [1.5:6]
set yrange [-3:10]
set grid 
replot

set output 'Conductivity.eps'
unset xrange
unset yrange
unset grid
set ylabel 'Conductivity {/Symbol s}'
p file u ($1*Hartree_to_eV):(real(conductivity($2,$3,Hartree_to_eV*$1/hbar*ec))) w l lt 1 lc "red" t 'Re({/Symbol s}_x)', \
  file u ($1*Hartree_to_eV):(imag(conductivity($2,$3,Hartree_to_eV*$1/hbar*ec))) w l lt 2 lc "blue" t 'Im({/Symbol s}_x)'

#
#  file u ($1*Hartree_to_eV):4 w l lt 1 lc 2 t 'Re({/Symbol e}_y)', \
#  file u ($1*Hartree_to_eV):5 w l lt 2 lc 2 t 'Im({/Symbol e}_y)', \
#  file u ($1*Hartree_to_eV):6 w l lt 1 lc 3 t 'Re({/Symbol e}_z)', \
#  file u ($1*Hartree_to_eV):7 w l lt 2 lc 3 t 'Im({/Symbol e}_z)'
