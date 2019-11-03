#!gnuplot

reset

set xrange [-10:25]

PathToFiles="/home/tjyderrien/octopus-slabs/06-Mo/"

H_to_eV=27.21138602 #1 Hartree -> eV

Zhigilei=PathToFiles."DOS_Mo.dat"
localDos="static/total-dos.dat"
Dorchies=PathToFiles."Dorchies2015-Fig1-Mo.csv"
Nechaev=PathToFiles."Nechaev-Fig1-Mo.csv"
Dos1="static/total-dos-3013.dat"
Dos2="static/total-dos-0418.dat"
Dos3="static/total-dos-4113.dat"
Dos4="static/total-dos-5183.dat"
EF1=3.012299 #Ha
EF2=0.418751 #Ha 
EF3=4.113481 #Ha
EF4=0.311046 #Ha
# localEF=0.462479 #Ha #PBE

## Finding latest LOG file automatically

localEFfile=system("ls -tr *.log -1 | tail -n1")
print "Latest LOG file: ".localEFfile
localEF=system("grep 'Fermi' ".localEFfile." | tail -n1 | awk '{ print $4 }'") #0.452100 #Ha with LDA. 
print "Fermi energy: ".localEF." H"

set xlabel 'Energy (eV)'

set key top left 

set ytics  nomirror
set y2tics nomirror
coeff=1./2.5/4./1.5

unset log y

p Zhigilei  u 1:2 w l t 'Lin and Zighilei (2007)' lc 1 lw 3, \
  localDos  u (($1-localEF)*H_to_eV):($2*coeff) w l t 'This work' axis x1y1 lt 3 lc 2 lw 3, \
  Nechaev u ($1):2 w l t 'Nechaev (2007)' axis x1y1 lc 4 lw 3

set terminal postscript eps enhanced color font 'Helvetica, 22'
set output 'DOS.eps'
replot

#  Dorchies u ($1-16.5):2 w l t 'Dorchies (2015)' axis x1y1 lc 3 lw 3, \


#  Dos2      u (($1-EF2)*H_to_eV):2 w l t 'ASE to Bohr' axis x1y2
#  Dos1      u (($1-EF1)*H_to_eV):2 w l t 'ASE (Angst.)' axis x1y2
#  Dos3      u (($1-EF3)*H_to_eV):2 w l t 'DOI: 10.17188/1189253 (Angst.)' axis x1y2
#  Dos4      u (($1-EF4)*H_to_eV):2 w l t 'DOI: 10.17188/1189253 (Bohr)' axis x1y2

#  localDos  u ($1*H_to_eV):2 w l t 'DOS' axis x1y2, \


