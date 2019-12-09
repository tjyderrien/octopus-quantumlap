#!gnuplot

reset
set terminal x11
set output

set xrange [-10:25]
Material="09-Au"
PathToFiles=system("echo $QuantumLaPsources")."/".Material."/"
# "/home/tjyderrien/octopus-slabs/09-Au/"

H_to_eV=27.21138602 #1 Hartree -> eV

Zhigilei=PathToFiles."DOS_Au.dat"
localDos="static/total-dos.dat"

## Finding latest LOG file automatically

localEFfile=system("ls -tr *.log -1 | tail -n1")
print "Latest LOG file: ".localEFfile
localEF=system("grep 'Fermi' ".localEFfile." | tail -n1 | awk '{ print $4 }'") #0.452100 #Ha with LDA. 
print "Fermi energy: ".localEF." H"

set xlabel 'Energy (eV)'

set key top left 

set ytics  nomirror
set y2tics nomirror
coeff=1. #/2.5/4./1.5

unset log y

p Zhigilei  u 1:2 w l t 'Lin and Zighilei (2007)' lc 1 lw 3, \
  localDos  u (($1-localEF)*H_to_eV):($2*coeff) w l t 'This work' axis x1y2 lt 3 lc 2 lw 3

set output 'DOS.eps'
set terminal postscript eps enhanced color font 'Helvetica, 22'

replot
