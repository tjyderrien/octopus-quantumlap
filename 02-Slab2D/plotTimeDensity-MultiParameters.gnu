#!gnuplot
reset
set output 'RoleOfPulseDuration-Density.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24' size 15cm,8cm

path="/ptmp/tjyderri/octopus-slabs.data/03-Bulk3D/Si-1x1x1-k-24x24x1-Mesh0.434637008855-Mesh0.434637008855-T0e0/"
filename='n_ex'

tau10_parameters="lambda800e-9-tau10e-15-field5E9-CEP0.0-StH2/lambda400e-9-tau20e-15-field0e0-CEP0.0-delay0e-15-StH/td.general/"
tau10_density=path.tau10_parameters.filename
print tau10_density
tau10_light=path.tau10_parameters.'laser'

tau20_parameters="lambda800e-9-tau20e-15-field5E9-CEP0.0-StH2/lambda400e-9-tau20e-15-field0e0-CEP0.0-delay0e-15-StH/td.general/"
tau20_density=path.tau20_parameters.filename

tau30_parameters="lambda800e-9-tau30e-15-field5E9-CEP0.0-StH2/lambda400e-9-tau20e-15-field0e0-CEP0.0-delay0e-15-StH/td.general/"
tau30_density=path.tau30_parameters.filename

set format "%g"

set xrange [:70] #tau50fs: 120, tau30fs: 70, tau10fs: 25

set key top left

set xlabel 'Time [fs]'
set ylabel 'Density [a.u]'
# set ylabel 'Vector field amplitude A(t) [a.u]'

# set xtics 2
set ytics nomirror
# set y2tics nomirror
set grid

xscale=2.74e-17*1E15

plot tau10_density u ($2*xscale):3 w l lw 1 t '{/Symbol t}=10 fs', \
     tau20_density u ($2*xscale):3 w l lw 1 t '{/Symbol t}=20 fs', \
     tau30_density u ($2*xscale):3 w l lw 1 t '{/Symbol t}=30 fs'

# "laser" u ($2*xscale):3 w l lc 1 lw 3 t 'Laser 1', \
# "laser" u ($2*xscale):6 w l lc 3 lw 3 t 'Laser 2'

# "total_density" u ($2*xscale):4 w p lw 1 lc 1 t 'Jy', \
# "total_density" u ($2*xscale):($5) w p lw 1 lc 3 t 'Jz', \
