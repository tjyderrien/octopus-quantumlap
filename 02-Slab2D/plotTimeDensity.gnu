#!gnuplot

reset

file='n_ex'

set output 'ExcitedElectrons.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'

set format "%g"

#set xrange [:25]

set xlabel 'Time [fs]'
set ylabel 'Density N_e [a.u]'

xscale=2.74e-17*1E15

plot file u ($2*xscale):3 w p lw 5 lc 7 t 'N_{cond}' 
#file u ($2*xscale):4 w l lw 5 lc 1 t 'Jy', \
#file u ($2*xscale):($5/10.) w l lw 5 lc 3 t 'Jz / 10'
