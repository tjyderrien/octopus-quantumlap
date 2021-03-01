#!gnuplot
reset
set output 'Energy.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24' 
#size 8,4

set key outside top center maxrow 3

set format "%g"

xmax=system("grep 'tauAsi =' ../inp | awk '{ print $3*1E15 }'")
print xmax

set xrange [:4*xmax] #tau50fs: 120, tau30fs: 70, tau10fs: 25
# set xrange [:0.5]

set ytics nomirror
set y2tics nomirror

set xlabel 'Time [fs]'
set ylabel 'Electron and Ion Energy [a.u.]'
# set y2label 'Electron Energy [a.u]'
set y2label 'Laser energy [a.u.]'

xscale=2.74e-17*1E15

plot \
"laser"  u ($2*xscale):3  w l lw 3 lc rgbcolor "red" t 'Laser' axis x1y2, \
"energy" u ($2*xscale):6  w l lw 5 lc rgbcolor "blue" t 'Electrons', \
"energy" u ($2*xscale):7  w l lw 5 lc rgbcolor "dark-yellow" lt 2 t 'Eigenvalues', \
"energy" u ($2*xscale):($8 == 0 ? NaN : $8) w lp lw 2 lc rgbcolor "dark-orange" t 'Hartree', \
"energy" u ($2*xscale):($9 == 0 ? NaN : $9) w lp lw 2 lc rgbcolor "dark-green" t 'int n(r) v_{xc}', \
"energy" u ($2*xscale):(($10+$11 == 0) ? NaN : ($10+$11)) w lp lw 2 lc rgbcolor "skyblue" t 'E_{xc}'
# "energy" u ($2*xscale):($7-$9-2*$8-$5) w lp lw 5 lc rgbcolor "red" t 'E_{kin}'
# "energy" u ($2*xscale):10 w p lw 2 lc rgbcolor "light-green" lt 2 t 'Exchange v_{x}', \
# "energy" u ($2*xscale):11 w p lw 2 lc rgbcolor "dark-yellow" lt 3 t 'Correlation v_{c}'

# "energy" u ($2*xscale):4 w l lw 5 lc 1 t 'Ions (kin)', \

# "energy" u ($2*xscale):3 w l lw 5 lc 7 t 'Total', \
