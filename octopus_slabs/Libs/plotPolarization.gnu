#!/bin/bash
reset

laserfile="laser"

set output "Polarization.eps"
set terminal postscript eps enhanced color font "Helvetica, 24" size 10cm, 10cm

# set xrange [-50:50]
extent=200
set yrange [-extent:extent]
set zrange [-extent:extent] 

# set view 70,20 #perspective
set view 90,90 #projection for all times
set size square # 4.,4.
set grid

# set xrange [0:1000]
unset xtics
set ytics 100
set ztics 100

splot \
laserfile u ($2):($3):($4) w l lc 1 lt 1 t 'Pulse 1', \
laserfile u ($2):($6):($7) w l lc 3 lt 1 t 'Pulse 2', \
laserfile u ($2):($3+$6):($4+$7) w l lc 7 lt 1 t 'Total'

