#!gnuplot

reset

# set xlabel 'Peak amplitude E_{peak} (V/m)'
# set ylabel 'Energy levels (eV)'

h=6.63e-34
c=3E8
e=1.6e-19

Egap=2.56

set xlabel 'V/nm'
#set ylabel '{/Symbol m}m'
set ylabel 'eV'

folder="/home/thibault/Documents/octopus/03-Bulk3D/Si-1x1x1-k-24x24x1-Mesh0.434637008855-Mesh0.434637008855-T0e0/"

hills="Peaks_"
valleys="Valleys_"
ext=".csv"

xunit=1E-9
yunit=1E-6

unset log x
unset log y
set xrange [0.4e9*xunit:6E9*xunit]
set yrange [:2.0]

set key outside above right
set xtics 1
set ytics 0.5

set terminal postscript eps enhanced color font 'Helvetica, 16' size 4cm, 4cm
set output "PositionOfPeaks-Hills.eps"

TDDFT_Hills    =folder."PositionOf".hills."WithErrorBars_TDDFT".ext
Keldysh_Hills  =folder."PositionOf".hills."WithErrorBars_Keldysh".ext
TDDFT_Valleys  =folder."PositionOf".valleys."TDDFT".ext
Keldysh_Valleys=folder."PositionOf".valleys."Keldysh".ext

set errorbars small

print "Plotting hills..."

p TDDFT_Hills   u 1:2:3  w errorbar t 'TDDFT'     lc "red"  lw 1 pt 1 ps 0.5, \
  Keldysh_Hills u 1:2:3  w errorbar t 'Keldysh'   lc "blue" lw 1 pt 2 ps 0.3

#   Keldysh_Hills u (xunit*$1):($3)       w p         notitle   lc "blue" lw 2 pt 2, \
#   Keldysh_Hills u (xunit*$1):($4)       w p         notitle   lc "blue" lw 2 pt 2, \
#   Keldysh_Hills u (xunit*$1):($5)       w p         notitle   lc "blue" lw 2 pt 2, \
#   Keldysh_Hills u (xunit*$1):($6)       w p         notitle   lc "blue" lw 2 pt 2, \
#   Keldysh_Hills u (xunit*$1):($7)       w p         notitle   lc "blue" lw 2 pt 2, \
#   Keldysh_Hills u (xunit*$1):($8)       w p         notitle   lc "blue" lw 2 pt 2, \
#   Keldysh_Hills u (xunit*$1):($9)       w p         notitle   lc "blue" lw 2 pt 2, \
#   Keldysh_Hills u (xunit*$1):($10)      w p         notitle   lc "blue" lw 2 pt 2
#   

unset errorbars

print "Plotting valleys..."

print TDDFT_Valleys

set xrange [0.4:3.5]
set yrange [:1.7]
set output "PositionOfPeaks-Valleys.eps"
set terminal postscript eps enhanced color font 'Helvetica, 16' size 4cm, 4cm

p TDDFT_Valleys   u ($1*xunit):2  w p t 'TDDFT'     lc "red"  lw 1 pt 1 ps 1.0, \
  Keldysh_Valleys u ($1*xunit):2  w p t 'Keldysh'   lc "blue" lw 1 pt 2 ps 1.0, \
  Egap/2.                w  l       notitle       lc "black" lw 0.5 lt 2, \
  Egap/3.                w  l       notitle       lc "black" lw 0.5 lt 2, \
  Egap/4.                w  l       notitle       lc "black" lw 0.5 lt 2
  




# p TDDFT_Valleys   u (xunit*$1):($2) w p t 'TDDFT' lc "red"    lw 9 lt 1, \
#   TDDFT_Valleys   u (xunit*$1):($3) w p           notitle lc "red"    lw 9 lt 1, \
#   TDDFT_Valleys   u (xunit*$1):($4) w p           notitle lc "red"    lw 9 lt 1, \
#   TDDFT_Valleys   u (xunit*$1):($5) w p           notitle lc "red"    lw 9 lt 1, \
#   TDDFT_Valleys   u (xunit*$1):($6) w p           notitle lc "red"    lw 9 lt 1, \
#   Keldysh_Valleys u (xunit*$1):($2) w p t 'Keldysh' lc "blue" lw 2 pt 2, \
#   Keldysh_Valleys u (xunit*$1):($3) w p           notitle lc "blue"   lw 2 pt 2, \
#   Keldysh_Valleys u (xunit*$1):($4) w p           notitle lc "blue"   lw 2 pt 2, \
#   Keldysh_Valleys u (xunit*$1):($5) w p           notitle lc "blue"   lw 2 pt 2, \
#   Keldysh_Valleys u (xunit*$1):($6) w p           notitle lc "blue"   lw 2 pt 2, \
#   Keldysh_Valleys u (xunit*$1):($7) w p           notitle lc "blue"   lw 2 pt 2
# 
