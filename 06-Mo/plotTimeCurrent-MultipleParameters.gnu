#!gnuplot
reset
# set output 'MultipleCurrent.eps'
# set terminal postscript eps enhanced color font 'Helvetica, 24' size 15cm,8cm

path="/scratch/work/user/tjyderrien/octopus-slabs.data/06-Mo/"
filename='total_current'

field="1E9"
tau20_parameters="Mo-1x1x1-k-24x24x1-Mesh0.434637008855-Mesh0.434637008855-T0e0/lambda800e-9-tau30E-15-field".field."-CEP0.0-SinSquared2/lambda400e-9-tau30E-15-field0E0-CEP0.0-delay0e-15-SinSquared/td.general/"
tau10_parameters="Mo-1x1x1-k-24x24x1-Mesh0.472431531364-Mesh0.472431531364-T0e0/lambda800e-9-tau30E-15-field".field."-CEP0.0-SinSquared2/lambda400e-9-tau30E-15-field0E0-CEP0.0-delay0e-15-SinSquared/td.general/"
tau30_parameters="Mo-1x1x1-k-24x24x1-Mesh0.283458918819-Mesh0.283458918819-T0e0/lambda800e-9-tau30E-15-field".field."-CEP0.0-SinSquared2/lambda400e-9-tau30E-15-field0E0-CEP0.0-delay0e-15-SinSquared/td.general/"

tau10_current=path.tau10_parameters.filename
tau10_light=path.tau20_parameters.'laser'
tau20_current=path.tau20_parameters.filename
tau30_current=path.tau30_parameters.filename

print tau10_light

set format "%g"

set xrange [30:40] #tau50fs: 120, tau30fs: 70, tau10fs: 25

set key top left

set xlabel 'Time [fs]'
set ylabel 'Current [a.u]'
set y2label 'Vector field amplitude A(t) [a.u]'

# set xtics 2
set ytics nomirror
set y2tics nomirror
set grid

xscale=2.74e-17*1E15

plot tau10_current u ($2*xscale):3 w l lw 1 t 'Mesh 0.48', \
     tau20_current u ($2*xscale):3 w l lw 1 t 'Mesh 0.43', \
     tau30_current u ($2*xscale):3 w l lw 1 t 'Mesh 0.28', \
     tau10_light u ($2*xscale):3 w l lc 1 lw 3 axis x1y2 t 'Laser'
# "laser" u ($2*xscale):6 w l lc 3 lw 3 t 'Laser 2'

# "total_current" u ($2*xscale):4 w p lw 1 lc 1 t 'Jy', \
# "total_current" u ($2*xscale):($5) w p lw 1 lc 3 t 'Jz', \
