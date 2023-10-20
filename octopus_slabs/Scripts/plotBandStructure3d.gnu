#!gnuplot

reset
print "We prepare band structure in 3D, then will cut it into a plane/line, to interpolate energy levels in that space."

# set pm3d
# set isosamples 10,10
coeff=3./8.
set xrange [coeff*0.9:coeff*1.1]
set yrange [coeff*0.9:coeff*1.1]


sp "< awk '{ if($3==0) print $1,$2,$4 }' Kozak-Si-TB09-k40-8unocc-BandStructure3D.dat | sort -g -k 1 | awk -f /home/thibault.derrien/Documents/flirt1d.fresh/Scripts/post/pm3d.awk" w p, \
   "< awk '{ if($3==0) print }' Kozak-Si-TB09-k40-8unocc-BandStructure3D.dat" u 1:2:5 w p, \
   "< awk '{ if($3==0) print }' Kozak-Si-TB09-k40-8unocc-BandStructure3D.dat" u 1:2:6 w p, \
   "< awk '{ if($3==0) print }' Kozak-Si-TB09-k40-8unocc-BandStructure3D.dat" u 1:2:7 w p, \
   "< awk '{ if($3==0) print }' Kozak-Si-TB09-k40-8unocc-BandStructure3D.dat" u 1:2:8 w p, \
   "< awk '{ if($3==0) print }' Kozak-Si-TB09-k40-8unocc-BandStructure3D.dat" u 1:2:9 w p, \
   "< awk '{ if($3==0) print }' Kozak-Si-TB09-k40-8unocc-BandStructure3D.dat" u 1:2:10 w p, \
   "< awk '{ if($3==0) print }' Kozak-Si-TB09-k40-8unocc-BandStructure3D.dat" u 1:2:11 w p, \
   "< awk '{ if($3==0) print }' Kozak-Si-TB09-k40-8unocc-BandStructure3D.dat" u 1:2:12 w p

