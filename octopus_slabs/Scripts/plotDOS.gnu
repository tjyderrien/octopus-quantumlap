#!gnuplot
reset
path="static/total-dos.dat"

set terminal postscript eps enhanced color font "Helvetica, 24"
set output 'DOS-EffectOfTe.eps'

eV=27.21138602
FermiLevel=system("tail -n1 Si-1x1x1-k-24x24x1-Angle0-Mesh0.4346370088553005-Mesh0.4346370088553005-T0e0-Te0.05eV-TB09/static/total-dos-efermi.dat | awk '{ print $1 }' | tr '\' ' '")
kb=1.38E-23; ec=1.6e-19
print "Fermi Level [eV] = "
print FermiLevel*eV

set xlabel 'Electron energy (eV)'
set ylabel 'DOS (arb. u.)'
set y2label 'Distribution'
set ytics nomirror
set y2tics nomirror

FermiDirac(x,y,Te)=1./(1.+exp((x-y)/(kb*Te/ec)))

p "Si-1x1x1-k-24x24x1-Angle0-Mesh0.4346370088553005-Mesh0.4346370088553005-T0e0-Te0.05eV-TB09/".path u ($1*eV):2 w l t '0.05 eV' ls 1 lc 1 axis x1y1, \
"Si-1x1x1-k-24x24x1-Angle0-Mesh0.4346370088553005-Mesh0.4346370088553005-T0e0-Te0.05eV-TB09/".path u ($1*eV):(FermiDirac($1*eV,FermiLevel*eV,0.05)) w l notitle ls 1 lc 1 axis x1y2, \
"Si-1x1x1-k-24x24x1-Angle0-Mesh0.4346370088553005-Mesh0.4346370088553005-T0e0-Te1.09eV-TB09/".path u ($1*eV):2 w l t '1.09 eV' ls 2 lc 2 axis x1y1, \
"Si-1x1x1-k-24x24x1-Angle0-Mesh0.4346370088553005-Mesh0.4346370088553005-T0e0-Te1.09eV-TB09/".path u ($1*eV):(FermiDirac($1*eV,FermiLevel*eV,1.09)) w l notitle ls 2 lc 2 axis x1y2, \
"Si-1x1x1-k-24x24x1-Angle0-Mesh0.4346370088553005-Mesh0.4346370088553005-T0e0-Te3.06eV-TB09/".path u ($1*eV):2 w l t '3.06 eV' ls 3 lc 3 axis x1y1, \
"Si-1x1x1-k-24x24x1-Angle0-Mesh0.4346370088553005-Mesh0.4346370088553005-T0e0-Te3.06eV-TB09/".path u ($1*eV):(FermiDirac($1*eV,FermiLevel*eV,3.06)) w l notitle ls 3 lc 3 axis x1y2



