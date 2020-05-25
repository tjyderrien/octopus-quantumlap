#!gnuplot
reset
set output 'Energy.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'

set format "%g"

xmax=system("grep 'tauAsi =' ../inp | awk '{ print $3*1E15 }'")
print xmax

set xrange [:3*xmax] #tau50fs: 120, tau30fs: 70, tau10fs: 25

set ytics nomirror
set y2tics nomirror

set xlabel 'Time [fs]'
set ylabel 'Ion Energy [a.u]'
set y2label 'Electron Energy [a.u]'

xscale=2.74e-17*1E15

plot "energy" u ($2*xscale):4 w l lw 5 lc 1 t 'Ions (kin)', \
"energy" u ($2*xscale):6 w l lw 5 lc 3 t 'Electrons' axis x1y2
# "energy" u ($2*xscale):3 w l lw 5 lc 7 t 'Total', \
