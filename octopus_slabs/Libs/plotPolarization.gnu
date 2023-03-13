#!/bin/bash
reset

laserfile="laser"

set output "Polarization.eps"
set terminal postscript eps enhanced color font "Helvetica, 24"


# set view 70,20 #perspective
set view 90,90 #projection for all times
set grid

set xrange [0:1000]

splot \
laserfile u ($2):($3):($4) w l lc 1 t 'Pulse 1', \
laserfile u ($2):($6):($7) w l lc 3 t 'Pulse 2', \
laserfile u ($2):($3+$6):($4+$7) w l lc 7 t 'Total'

