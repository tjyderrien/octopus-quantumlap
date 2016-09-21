#!gnuplot
reset
set output 'Laser.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'

set format "%g"

#set xrange [:25]

set xlabel 'Time [fs]'
set ylabel 'Energy [a.u]'

xscale=2.74e-17*1E15

plot "laser" u ($2*xscale):3 w l lw 5 lc 1 t 'Laser 1', \
"laser" u ($2*xscale):6 w l lw 5 lc 3 t 'Laser 2'
