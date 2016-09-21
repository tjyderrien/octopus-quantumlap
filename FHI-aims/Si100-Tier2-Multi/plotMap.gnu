#!/bin/bash
reset

print "Selecting Samples with dimX=dimY=1..."
print "30 A vacuum give same results than 50 A vacuum. Merging data..." 

Selection1="< awk '{ if($2==1 && $4==30) print }' energy-vs-kgrid-filtered.dat | awk -f '../../Video/pm3d.awk'"
Selection2="< awk '{ if($2==1 && $4==50) print }' energy-vs-kgrid-filtered.dat | awk -f '../../Video/pm3d.awk'"

Selection="< awk '{ if($2==1 && $3==8) print }' energy-vs-kgrid-filtered.dat"
# PLOT 1
print "Plot 3d..."
set terminal postscript eps enhanced color font 'Times, 22'
set output 'Convergence.eps'
load '../../Video/RdYlBu.plt'

set xlabel 'Layers'
set ylabel 'k-points'
set zlabel 'Energy (eV)'

set xrange [4:16]
set yrange [4:16]
# set ztics 0.1
# set ctics

set xtics 2
set ytics 2
# set cbrange
set log xy
# set size 0.9
# set key autotitle columnheader
# set pm3d interpolate 8,8
set view 60,30

set cbtics 0.1
# unset cbtics
# set view map
splot Selection u 1:3:5 w p notitle

############################
# PLOT 2
#########
reset
print "Effect of vacuum size... [E(layers; vacuum)"

set output 'EffectOfVacuumLayers.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'
set xlabel 'Layers'
set ylabel 'Energy'

set xrange [:]
set yrange [:]
set xtics auto
set ytics auto
plot Selection1 u 1:5 w p t 'Vacuum: 30', \
Selection2 u 1:5 w p lc 3 t 'Vacuum: 50'

#### Plot 3 ###
reset

set output 'NumberOfLayers.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'
print "/!\\ 50 and 30 A vacuum are assumed equal."
print "/!\\ Assuming that 8 k-points are enough..."

set xlabel 'Number of Si cells in depth'
set ylabel 'Total energy per atom (eV)'

plot Selection u 1:5 w lp lw 3 lc 3 t 'Energy'
