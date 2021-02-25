#!gnuplot

reset

# set xlabel 'Peak amplitude E_{peak} (V/m)'
# set ylabel 'Energy levels (eV)'

h=6.63e-34
c=3E8
e=1.6e-19

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
# set xrange [0.4e9*xunit:4E9*xunit]
# set yrange [:1.2]


set key outside above right
set xtics 1
set ytics 0.5

# set terminal postscript eps enhanced color font 'Helvetica, 16' size 4cm, 4cm
# set output "PositionOfPeaks-Hills.eps"

TDDFT_Hills    =folder."PositionOf".hills."TDDFT".ext
Keldysh_Hills  =folder."PositionOf".hills."Keldysh".ext
TDDFT_Valleys  =folder."PositionOf".valleys."TDDFT".ext
Keldysh_Valleys=folder."PositionOf".valleys."Keldysh".ext

print "Plotting hills..."

p TDDFT_Hills   u (xunit*$1):($2)  w p t 'TDDFT'   lc "red"  lw 9 lt 1, \
  TDDFT_Hills   u (xunit*$1):($3)  w p         notitle   lc "red"  lw 9 lt 1, \
  TDDFT_Hills   u (xunit*$1):($4)  w p         notitle   lc "red"  lw 9 lt 1, \
  TDDFT_Hills   u (xunit*$1):($5)  w p         notitle   lc "red"  lw 9 lt 1, \
  TDDFT_Hills   u (xunit*$1):($6)  w p         notitle   lc "red"  lw 9 lt 1, \
  TDDFT_Hills   u (xunit*$1):($7)  w p         notitle   lc "red"  lw 9 lt 1, \
  TDDFT_Hills   u (xunit*$1):($8)  w p         notitle   lc "red"  lw 9 lt 1, \
  Keldysh_Hills u (xunit*$1):($2)  w p t 'Keldysh' lc "blue" lw 2 pt 2, \
  Keldysh_Hills u (xunit*$1):($3)  w p         notitle   lc "blue" lw 2 pt 2, \
  Keldysh_Hills u (xunit*$1):($4)  w p         notitle   lc "blue" lw 2 pt 2, \
  Keldysh_Hills u (xunit*$1):($5)  w p         notitle   lc "blue" lw 2 pt 2, \
  Keldysh_Hills u (xunit*$1):($6)  w p         notitle   lc "blue" lw 2 pt 2, \
  Keldysh_Hills u (xunit*$1):($7)  w p         notitle   lc "blue" lw 2 pt 2, \
  Keldysh_Hills u (xunit*$1):($8)  w p         notitle   lc "blue" lw 2 pt 2, \
  Keldysh_Hills u (xunit*$1):($9)  w p         notitle   lc "blue" lw 2 pt 2, \
  Keldysh_Hills u (xunit*$1):($10) w p         notitle   lc "blue" lw 2 pt 2
  
print "Plotting valleys..."

set terminal x11 2
# set output "PositionOfPeaks-Valleys.eps"

p TDDFT_Valleys   u (xunit*$1):($2) w p t 'TDDFT' lc "red"    lw 9 lt 1, \
  TDDFT_Valleys   u (xunit*$1):($3) w p           notitle lc "red"    lw 9 lt 1, \
  TDDFT_Valleys   u (xunit*$1):($4) w p           notitle lc "red"    lw 9 lt 1, \
  TDDFT_Valleys   u (xunit*$1):($5) w p           notitle lc "red"    lw 9 lt 1, \
  TDDFT_Valleys   u (xunit*$1):($6) w p           notitle lc "red"    lw 9 lt 1, \
  Keldysh_Valleys u (xunit*$1):($2) w p t 'Keldysh' lc "blue" lw 2 pt 2, \
  Keldysh_Valleys u (xunit*$1):($3) w p           notitle lc "blue"   lw 2 pt 2, \
  Keldysh_Valleys u (xunit*$1):($4) w p           notitle lc "blue"   lw 2 pt 2, \
  Keldysh_Valleys u (xunit*$1):($5) w p           notitle lc "blue"   lw 2 pt 2, \
  Keldysh_Valleys u (xunit*$1):($6) w p           notitle lc "blue"   lw 2 pt 2, \
  Keldysh_Valleys u (xunit*$1):($7) w p           notitle lc "blue"   lw 2 pt 2
