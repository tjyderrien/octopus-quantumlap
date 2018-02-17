#!gnuplot

reset

nex='n_ex'
light='laser'


set output 'ExcitedElectrons.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'

set format "%g"

set xrange [:70]

set xlabel 'Time [fs]'
set ylabel 'A(t) (a.u.)'
set y2label 'Density N_e [a.u]'

set ytics nomirror
set y2tics nomirror

set key right bottom

xscale=2.74e-17*1E15

plot nex u ($2*xscale):3 w l lw 5 lc 0 axis x1y2 t 'N_{cond}', \
light u ($2*xscale):(($3**2)**0.5) w l lw 1 lc 1 lt 2 t 'A_1(t)', \
light u ($2*xscale):(($6**2)**0.5) w l lw 1 lc 2 lt 3 t 'A_2(t)', \
light u ($2*xscale):((($3+$6)**2)**0.5) w l lw 5 lc 3 t 'A_{tot}(t)'
