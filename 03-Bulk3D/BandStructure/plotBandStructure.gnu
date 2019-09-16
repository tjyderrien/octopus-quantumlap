#!gnuplot

# set terminal x11

# set output 'BandStructure.eps'
# set terminal postscript eps enhanced font 'Helvetica, 24'
file="bandstructure"

HartreeToEv = 27.2114 #1 Hartree = 27.2114 eV

Yscale = HartreeToEv

#Column meanings. 
# 1 2  3  4  5     6     7     ... 8 9 
# ? kx ky kz band1 band2 band3 ....

plot file u 1:($5*Yscale) w  l notitle, \
     file u 1:($6*Yscale) w  l notitle, \
     file u 1:($7*Yscale) w  l notitle, \
     file u 1:($8*Yscale) w  l notitle, \
     file u 1:($9*Yscale) w  l notitle, \
     file u 1:($10*Yscale) w l notitle, \
     file u 1:($11*Yscale) w l notitle, \
     file u 1:($12*Yscale) w l notitle, \
     file u 1:($13*Yscale) w l notitle, \
     file u 1:($14*Yscale) w l notitle




