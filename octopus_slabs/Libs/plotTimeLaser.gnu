#!gnuplot
reset
set output 'Laser.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'

set format "%g"

# set xrange [:50]

set ytics nomirror
set y2tics nomirror

set xlabel 'Time t [fs]'
set ylabel 'Vector field amplitude A_1(t) [a.u]'
set y2label 'Vector field ampliude A_2(t) [a.u]'

xscale=2.74e-17*1E15

plot "laser.tmp" u ($2*xscale):3 w l lw 5 lc 1 t 'A1' axis x1y1, \
     "laser.tmp" u ($2*xscale):6 w l lw 5 lc 3 t 'A2' axis x1y2

