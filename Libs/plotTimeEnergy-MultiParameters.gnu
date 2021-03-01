#!gnuplot
reset

path="/ptmp/tjyderri/octopus-slabs.data/03-Bulk3D/Si-1x1x1-k-24x24x1-Mesh0.434637008855-Mesh0.434637008855-T0e0/"

filename='energy'

# field="1E10" #1E10
pulse10="10e-15"
pulse20="20e-15"
pulse30="30E-15"

set output 'RoleOfPulseDuration-Energy'.field.'.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24' size 15cm,8cm

tau10_parameters="lambda800e-9-tau".pulse10."-field".field."-CEP0.0-StH2/lambda400e-9-tau20e-15-field0e0-CEP0.0-delay0e-15-StH/td.general/"
tau10_energy=path.tau10_parameters.filename
print tau10_energy
tau10_light=path.tau10_parameters.'laser'

tau20_parameters="lambda800e-9-tau".pulse20."-field".field."-CEP0.0-StH2/lambda400e-9-tau20e-15-field0e0-CEP0.0-delay0e-15-StH/td.general/"
tau20_energy=path.tau20_parameters.filename

tau30_parameters="lambda800e-9-tau".pulse30."-field".field."-CEP0.0-StH2/lambda400e-9-tau20E-15-field0E0-CEP0.0-delay0e-15-StH/td.general/"
tau30_energy=path.tau30_parameters.filename

set format "%g"

set xrange [:50] #tau50fs: 120, tau30fs: 70, tau10fs: 25
set yrange [0:] 
set key top left

set xlabel 'Time [fs]'
unset ylabel #set ylabel 'Energy (Ha per atom)'
# set ylabel 'Vector field amplitude A(t) [a.u]'

# set xtics 2
set ytics nomirror
# set y2tics nomirror
unset grid

NumberOfAtoms = 2

xscale=2.74e-17*1E15
Origin=-0.06+0.00092

plot tau10_energy u ($2*xscale):($6/NumberOfAtoms-Origin) w l lc rgb "black" lw 3 t '{/Symbol t}=10 fs', \
     tau20_energy u ($2*xscale):($6/NumberOfAtoms-Origin) w l lc rgb "blue"  lw 3 t '{/Symbol t}=20 fs', \
     tau30_energy u ($2*xscale):($6/NumberOfAtoms-Origin) w l lc rgb "red"   lw 3 t '{/Symbol t}=30 fs'

# "laser" u ($2*xscale):3 w l lc 1 lw 3 t 'Laser 1', \
# "laser" u ($2*xscale):6 w l lc 3 lw 3 t 'Laser 2'

# "total_energy" u ($2*xscale):4 w p lw 1 lc 1 t 'Jy', \
# "total_energy" u ($2*xscale):($5) w p lw 1 lc 3 t 'Jz', \
