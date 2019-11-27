#!gnuplot

reset

# set xlabel 'Peak amplitude E_{peak} (V/m)'
# set ylabel 'Energy levels (eV)'

set xlabel 'V/m'
set ylabel '{/Symbol m}m'

folder="/home/thibault/Documents/octopus/03-Bulk3D/Si-1x1x1-k-24x24x1-Mesh0.434637008855-Mesh0.434637008855-T0e0/"

hills="Peaks_"
valleys="Valleys_"
ext=".csv"

set log x
set xrange [0.4e9:4E9]
set yrange [:1.2]

set key outside below center

set terminal postscript eps enhanced color font 'Helvetica, 16' size 4cm, 4cm

set output "PositionOfPeaks-Hills.eps"

TDDFT_Hills    =folder."PositionOf".hills."TDDFT".ext
Keldysh_Hills  =folder."PositionOf".hills."Keldysh".ext
TDDFT_Valleys  =folder."PositionOf".valleys."TDDFT".ext
Keldysh_Valleys=folder."PositionOf".valleys."Keldysh".ext

print "Plotting hills..."

p TDDFT_Hills u 1:2 w p t 'TDDFT Hills' lc "red"    lw 9 lt 1, \
TDDFT_Hills u 1:3 w p           notitle lc "red"    lw 9 lt 1, \
TDDFT_Hills u 1:4 w p           notitle lc "red"    lw 9 lt 1, \
TDDFT_Hills u 1:5 w p           notitle lc "red"    lw 9 lt 1, \
TDDFT_Hills u 1:6 w p           notitle lc "red"    lw 9 lt 1, \
Keldysh_Hills u 1:2 w p t 'Keldysh Hills' lc "blue" lw 5 pt 2, \
Keldysh_Hills u 1:3 w p           notitle lc "blue" lw 5 pt 2, \
Keldysh_Hills u 1:4 w p           notitle lc "blue" lw 5 pt 2, \
Keldysh_Hills u 1:5 w p           notitle lc "blue" lw 5 pt 2, \
Keldysh_Hills u 1:6 w p           notitle lc "blue" lw 5 pt 2
  
print "Plotting valleys..."

# set terminal x11 2
set output "PositionOfPeaks-Valleys.eps"

p TDDFT_Valleys u 1:2 w p t 'TDDFT Valleys' lc "red"      lw 9 lt 1, \
  TDDFT_Valleys u 1:3 w p           notitle lc "red"      lw 9 lt 1, \
  TDDFT_Valleys u 1:4 w p           notitle lc "red"      lw 9 lt 1, \
  TDDFT_Valleys u 1:5 w p           notitle lc "red"      lw 9 lt 1, \
  TDDFT_Valleys u 1:6 w p           notitle lc "red"      lw 9 lt 1, \
  Keldysh_Valleys u 1:2 w p t 'Keldysh Valleys' lc "blue" lw 5 pt 2, \
  Keldysh_Valleys u 1:3 w p             notitle lc "blue" lw 5 pt 2, \
  Keldysh_Valleys u 1:4 w p             notitle lc "blue" lw 5 pt 2, \
  Keldysh_Valleys u 1:5 w p             notitle lc "blue" lw 5 pt 2, \
  Keldysh_Valleys u 1:6 w p             notitle lc "blue" lw 5 pt 2
