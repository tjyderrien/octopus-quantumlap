#!/bin/bash
reset

print "Selecting Samples with dimX=dimY=1..."
print "30 A vacuum give same results than 50 A vacuum. Merging data..." 

DataFile="energy-vs-kgrid-filtered.dat"

Selection1="< awk '{ if($2==1 && $4==30) print }' ".DataFile."| awk -f '../../Video/pm3d.awk'"
Selection2="< awk '{ if($2==1 && $4==50) print }' ".DataFile." | awk -f '../../Video/pm3d.awk'"
Selection="< awk '{ if($2==1 && $3==8) print }' ".DataFile

Slab1="< awk '{ if($3==8) print }' ".DataFile." | awk -f '../../Video/pm3d.awk'"
Slab2="< awk '{ if($3==1) print }' ".DataFile." | awk -f '../../Video/pm3d.awk'"

IncreasingZ="< awk '{ if($2==1 && $3==8) print }' ".DataFile
IncreasingX16="< awk '{ if($1==10 && $3==16) print }' ".DataFile
IncreasingX8="< awk '{ if($1==10 && $3==8) print }' ".DataFile
IncreasingX1="< awk '{ if($1==10 && $3==1) print }' ".DataFile
# Slab21="< awk '{ if($2==2 && $3==1) print }' DataFile"

# Tier2="< awk '{ if($2==1 && $3==8) print }' ../Si100-Tier2-Multi/DataFile"

load '../../Video/RdYlBu.plt'

############### PLOT 3D: Energy (k-points, layersZ)
set terminal postscript eps enhanced color font 'Times, 22'
set output 'Convergence.eps'
print "Plot 3d..."
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

#########
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

#### Plot 3 ####
# reset
# 
# set output 'NumberOfLayers.eps'
# set terminal postscript eps enhanced color font 'Helvetica, 24'
# print "/!\\ 50 and 30 A vacuum are assumed equal."
# print "/!\\ Assuming that 8 k-points are enough..."
# 
# set xlabel 'Number of Si cells in depth'
# set ylabel 'Total energy per atom (eV)'
# 
# plot Selection u 1:5 w lp lw 3 lc 1 t 'Tier1', \
# Tier2 u 1:5 w lp lw 3 lc 3 t 'Tier2'

############### PLOT4: Energy(LayerZ) for two cases. 
# 1. LayerX = 1 and k=1x8x8
# 2. Layer X increases and k=1x1x1. 
print "Plotting Energy(LayerZ, LayerXY)"
reset

set output 'Convergence-in-space.eps'
set terminal postscript eps enhanced color font 'Times, 22'

set xlabel 'N_{depth}'
set ylabel 'N_{||}'

set title 'Total energy per atom (eV)'

# plot Slab1 u 1:5 w l lp lw 3 lc 1 t 'Slab1: k=1x8x8',
# Slab21 u 1:5 w l lp lw 3 lc 3 t 'Slab2: XxY=1x1, k=1x1x1'
Eshift=-7159.75
# set view map
splot Slab1 u 1:2:($5-Eshift) w p lc 1 t 'Slab1: k=1x8x8', \
Slab2 u 1:2:($5-Eshift) w p lc 3 t 'Slab2: k=1x1x1' 

# Slab22 u 1:5 w l lp lw 3 lc 3 t 'Slab2: XxY=2x2, k=1x1x1',
# Slab23 u 1:5 w l lp lw 3 lc 3 t 'Slab2: XxY=3x3, k=1x1x1', 
# Slab24 u 1:5 w l lp lw 3 lc 3 t 'Slab2: XxY=4x4, k=1x1x1'

###### PLOT5: let's make it simpler. 
# 1. Energy(LayerZ; k=1x8x8)
# 2. Energy(LayerX; k=1x8x8) and Energy(LayerX; k=1x1x1). 

print "Plotting Energy(LayerZ; k=1x8x8)..."

reset

set output 'Convergence-LayerZ.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'
# set terminal x11
# set output 

set xlabel 'Number of blocks in depth'
set ylabel 'Surface energy per cell (eV)'

plot IncreasingZ u 1:7 w lp lw 3 t 'k=1x8x8'

print "Plotting Energy(LayerX)..."
set output 'Convergence-LayerX.eps'
set xlabel 'Number of cells along the surface'
set ylabel 'Surface energy per cell (eV)'

plot IncreasingZ u 2:7 w lp lw 3 lc 1 t 'k=1x8x8' #, \
# IncreasingX1 u 2:7 w lp lw 3 lc 3 t 'k=1x1x1'
