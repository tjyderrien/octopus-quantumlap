#!gnuplot
reset
set output 'Current.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24' size 15cm,8cm

set format "%g"

set xrange [:25] #tau50fs: 120, tau30fs: 70, tau10fs: 25

set key top left

set xlabel 'Time [fs]'
set ylabel 'Vector field amplitude A(t) [a.u]'
set y2label 'Current [a.u]'

# set xtics 2
set ytics nomirror
set y2tics nomirror
set grid

xscale=2.74e-17*1E15

plot "total_current" u ($2*xscale):3 w l lw 5 lc 7 t 'Jx' axis x1y2, \
"laser" u ($2*xscale):3 w l lc 1 lw 3 t 'Laser 1', \
"laser" u ($2*xscale):6 w l lc 3 lw 3 t 'Laser 2'

# "total_current" u ($2*xscale):4 w p lw 1 lc 1 t 'Jy', \
# "total_current" u ($2*xscale):($5) w p lw 1 lc 3 t 'Jz', \
