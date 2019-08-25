#!gnuplot

reset

set xrange [-10:25]

H_to_eV=27.21138602 #1 Hartree -> eV

Zhigilei="DOS_Mo.dat"
localDos="static/total-dos.dat"
Dorchies="Dorchies2015-Fig1-Mo.csv"
Nechaev="Nechaev-Fig1-Mo.csv"
Dos1="static/total-dos-3013.dat"
Dos2="static/total-dos-0418.dat"
Dos3="static/total-dos-4113.dat"
Dos4="static/total-dos-5183.dat"
EF1=3.012299 #Ha
EF2=0.418751 #Ha 
EF3=4.113481 #Ha
EF4=0.311046 #Ha
localEF=0.462062 #Ha
set xlabel 'Energy (eV)'

set ytics  nomirror
set y2tics nomirror

# set y2range [:5]

p Zhigilei  u 1:2 w l t 'Zighilei', \
  localDos  u (($1-localEF)*H_to_eV):2 w lp t 'local' axis x1y2, \
  Dorchies u ($1-16.5):2 w l t 'Dorchies (2015)' axis x1y1, \
  Nechaev u ($1):2 w l t 'Nechaev (2007)' axis x1y1

#  Dos2      u (($1-EF2)*H_to_eV):2 w l t 'ASE to Bohr' axis x1y2
#  Dos1      u (($1-EF1)*H_to_eV):2 w l t 'ASE (Angst.)' axis x1y2
#  Dos3      u (($1-EF3)*H_to_eV):2 w l t 'DOI: 10.17188/1189253 (Angst.)' axis x1y2
#  Dos4      u (($1-EF4)*H_to_eV):2 w l t 'DOI: 10.17188/1189253 (Bohr)' axis x1y2




#  localDos  u ($1*H_to_eV):2 w l t 'DOS' axis x1y2, \


