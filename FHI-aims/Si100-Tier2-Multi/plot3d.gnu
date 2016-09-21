#!gnuplot

reset

# set terminal postscript eps enhanced color font 'Times, 22'
# set output 'Convergence.eps'
load '../../Video/RdYlBu.plt'

set xlabel 'Layers'
set ylabel 'k-points'
set zlabel 'Energy (eV)'

# set xrange [4:16]
# set yrange [4:8]
# set ctics

set xtics 2
set ytics 2
# set cbrange
set log xy
# set size 0.9
# set key autotitle columnheader
# set pm3d interpolate 8,8

# set view map
splot "< awk -f '../../Video/pm3d.awk' energy-vs-kgrid-filtered.dat" u 1:2:4 w p notitle

