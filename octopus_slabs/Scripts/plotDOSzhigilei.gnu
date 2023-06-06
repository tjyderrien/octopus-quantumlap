#!gnuplot

reset
set terminal x11
set output

set xrange [-10:25]

Material="Au" #Mo|Au
#PathToFiles="/home/tjyderrien/octopus_slabs/06-Mo/"
PathToFiles="/home/tjyderrien/octopus_slabs/09-Au/"

H_to_eV=27.21138602 #1 Hartree -> eV

Zhigilei=PathToFiles."DOS_".Material.".dat"
localDos="static/total-dos.dat"
Dorchies=PathToFiles."Dorchies2015-Fig1-Mo.csv"
Nechaev=PathToFiles."Nechaev-Fig1-Mo.csv"

## Finding latest LOG file automatically

localEFfile=system("ls -tr *octopus*.log -1 | tail -n1")
print "Latest LOG file: ".localEFfile
localEF=system("grep 'Fermi' ".localEFfile." | tail -n1 | awk '{ print $4 }'") #0.452100 #Ha with LDA. 
print "Fermi energy: ".localEF." H"

set xlabel 'Energy (eV)'

set key top left 

# set ytics  nomirror
# set y2tics nomirror
coeff=1./2.5/4./1.5

unset log y

p Zhigilei  u 1:2 w l t 'Lin and Zighilei (2007)' lc 1 lw 3, \
  localDos  u (($1-localEF)*H_to_eV):($2*coeff) w l t 'This work' axis x1y1 lt 3 lc 2 lw 3, \
  Nechaev u ($1):2 w l t 'Nechaev (2007)' axis x1y1 lc 4 lw 3

set output 'DOS.eps'
set terminal postscript eps enhanced color font 'Helvetica, 22'

replot
