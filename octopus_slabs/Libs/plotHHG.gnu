#!gnuplot

reset

field="3E9"

set terminal postscript eps enhanced color font "Helvetica, 20"
set output 'Kozak_HHG_Si_'.field.'Vm.eps'


set xlabel 'Energy (eV)'
set ylabel 'Intensity (a. u.)'

set log y

set xrange [0.1:10]

path1="Si-1x1x1-k-24x24x1-Angle"
path2="-Mesh0.434637008855-Mesh0.434637008855-T0e0/lambda2656e-9-tau30e-15-field"
path3="-CEP0.0-theta"
path4="-SinSquared-SinPower2-2/lambda400e-9-tau30e-15-field0E0-CEP0.0-delay0e-15-SinSquared/td.general/"
path5="HHG-Spectrum.csv"

p path1."0".path2.field.path3."0".path4.path5 w l t "0", \
path1."0.05".path2.field.path3."0.05".path4.path5 w l t "0.05", \
path1."0.10".path2.field.path3."0.10".path4.path5 w l t "0.10", \
path1."0.15".path2.field.path3."0.15".path4.path5 w l t "0.15", \
path1."0.20".path2.field.path3."0.20".path4.path5 w l t "0.20", \
path1."0.25".path2.field.path3."0.25".path4.path5 w l t "0.25"

###==========================================================

set output 'Kozak_Currents_'.field.'Vm.eps'
set xlabel 'Time (at. u.)'
set ylabel 'Current (at. u.)'

unset log xy
set xrange [0:2800]

p path1."0".path2.field.path3."0".path4."total_current"     u 2:3 every 10 w l  lc 1 lt 1 t "Jx: 0", \
path1."0.05".path2.field.path3."0.05".path4."total_current" u 2:3 every 10 w l  lc 2 lt 1 t "Jx: 0.05 pi", \
path1."0.10".path2.field.path3."0.10".path4."total_current" u 2:3 every 10 w l  lc 3 lt 1 t "Jx: 0.10 pi", \
path1."0.15".path2.field.path3."0.15".path4."total_current" u 2:3 every 10 w l  lc 4 lt 1 t "Jx: 0.15 pi", \
path1."0.20".path2.field.path3."0.20".path4."total_current" u 2:3 every 10 w l  lc 5 lt 1 t "Jx: 0.20 pi", \
path1."0.25".path2.field.path3."0.25".path4."total_current" u 2:3 every 10 w l  lc 6 lt 1 t "Jx: 0.25 pi", \
path1."0".path2.field.path3."0".path4."total_current"       u 2:4 every 10 w lp lc 1 lt 2 t "Jy: 0", \
path1."0.05".path2.field.path3."0.05".path4."total_current" u 2:4 every 10 w lp lc 2 lt 2 t "Jy: 0.05 pi", \
path1."0.10".path2.field.path3."0.10".path4."total_current" u 2:4 every 10 w lp lc 3 lt 2 t "Jy: 0.10 pi", \
path1."0.15".path2.field.path3."0.15".path4."total_current" u 2:4 every 10 w lp lc 4 lt 2 t "Jy: 0.15 pi", \
path1."0.20".path2.field.path3."0.20".path4."total_current" u 2:4 every 10 w lp lc 5 lt 2 t "Jy: 0.20 pi", \
path1."0.25".path2.field.path3."0.25".path4."total_current" u 2:4 every 10 w lp lc 6 lt 2 t "Jy: 0.25 pi"
