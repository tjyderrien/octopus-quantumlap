#!gnuplot
reset

set output 'Summary-ElectronDensityWithIntensity.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'

set xlabel 'Laser intensity (W/cm^{2})'
set ylabel 'Ionization probability w_{PI} (m^{-3} s^{-1})'

DatabaseK24="/ptmp/tjyderri/octopus-slabs.data/03-Bulk3D/Si-1x1x1-k-24x24x1-Mesh0.434637008855-Mesh0.434637008855-T0e0/Database.csv"
DatabaseK16="/ptmp/tjyderri/octopus-slabs.data/03-Bulk3D/Si-1x1x1-k-16x16x1-Mesh0.434637008855-Mesh0.434637008855-T0e0/Database.csv"



