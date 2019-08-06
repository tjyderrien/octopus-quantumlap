#!gnuplot

reset

H_to_eV=27.21138602 #1 Hartree -> eV

Zhigilei="DOS_Mo.dat"
localDos="static/total-dos.dat"
Dorchies="Dorchies2015-Fig1-Mo.csv"
Dos1="static/total-dos-3147.dat"
Dos2="static/total-dos-5946.dat"
Dos3="static/total-dos-2743.dat"
Dos4="static/total-dos-5183.dat"
set xlabel 'Energy (eV)'

set ytics  nomirror
set y2tics nomirror

# set y2range [:5]

p Zhigilei  u 1:2 w l t 'Zighilei', \
  Dos4      u ($1*H_to_eV):2 w l t 'DOI: 10.17188/1189253 (Bohr)' axis x1y2
#  Dos1      u ($1*H_to_eV):2 w l t 'ASE (Angst.)' axis x1y2
#  Dos2      u ($1*H_to_eV):2 w l t 'ASE to Bohr' axis x1y2
#  Dos3      u ($1*H_to_eV):2 w l t 'DOI: 10.17188/1189253 (Angst.)' axis x1y2, \


#  localDos  u ($1*H_to_eV):2 w l t 'DOS' axis x1y2, \

#  Dorchies u ($1-16.5):2 w l t 'Dorchies (2015)' axis x1y1
