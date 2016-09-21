#!gnuplot
reset
set output 'Current.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'

set format "%g"

#set xrange [:2]

set xlabel 'Time [fs]'
set ylabel 'Current [a.u]'

xscale=2.74e-17*1E15

plot "total_current" u ($2*xscale):3 w p lw 3 lc 7 t 'Jx', \
"total_current" u ($2*xscale):4 w p lw 3 lc 1 t 'Jy', \
"total_current" u ($2*xscale):($5) w p lw 3 lc 3 t 'Jz'
