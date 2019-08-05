#!gnuplot

reset

H_to_eV=27.21138602 #1 Hartree -> eV

base="DOS_Mo.dat"
localDos="static/total-dos.dat"

set xlabel 'Energy (eV)'

p base  u 1:2 w l t 'Zighilei', \
  localDos u ($1*H_to_eV):2 w l t 'DOS' axis x1y2
