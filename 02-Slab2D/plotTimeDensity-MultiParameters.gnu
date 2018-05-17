#!gnuplot

reset
path="/ptmp/tjyderri/octopus-slabs.data/03-Bulk3D/Si-1x1x1-k-16x16x1-Mesh0.434637008855-Mesh0.434637008855-T0e0/"
tau10_parameters="lambda800e-9-tau10e-15-field5e9-CEP0.0-STH/lambda400e-9-tau10e-15-field0e9-CEP0.0-delay0e-15-STH/td.general/"
tau10_nex=path.tau10_parameters.'n_ex'
tau10_light=path.tau10_parameters.'laser'

tau20_parameters="lambda800e-9-tau20e-15-field5e9-CEP0.0-STH/lambda400e-9-tau10e-15-field0e9-CEP0.0-delay0e-15-STH/td.general/"
tau20_nex=path.tau20_parameters.'n_ex'

tau30_parameters="lambda800e-9-tau30e-15-field5e9-CEP0.0-STH/lambda400e-9-tau10e-15-field0e9-CEP0.0-delay0e-15-STH/td.general/"
tau30_nex=path.tau30_parameters.'n_ex'

set output 'RoleOfPulseDuration.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'

set format "%g"

set xrange [-5:50]

set xlabel 'Time [fs]'
# set y2label 'A(t) (a.u.)'
set ylabel 'Density N_e [a.u]'

set ytics nomirror
# set y2tics nomirror

set key right bottom

xscale=2.74e-17*1E15

## We compare different pulse durations for the same energy. 

plot tau10_nex u ($2*xscale-5):3 w l lw 3 lc 1 t '{/Symbol t}=10 fs', \
tau20_nex u ($2*xscale-10):3 w l lw 3 lc 2 t '{/Symbol t}=20 fs', \
tau30_nex u ($2*xscale-15):3 w l lw 3 lc 3 t '{/SYmbol t}=30 fs'

# light u ($2*xscale):(($3**2)**0.5) w l lw 1 lc 1 lt 2 t 'A(t)'
                         
