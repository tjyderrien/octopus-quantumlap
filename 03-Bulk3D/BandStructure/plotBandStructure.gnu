#!gnuplot

# set terminal x11

set terminal postscript eps enhanced font 'Helvetica, 24'
set output 'BandStructure.eps'


file="bandstructure"

HartreeToEv = 27.2114 #1 Hartree = 27.2114 eV

Yscale = HartreeToEv

set ylabel 'Energy (eV)'

unset xtics

#Column meanings. 
# 1 2  3  4  5     6     7     ... 8 9 
# ? kx ky kz band1 band2 band3 ....


FermiLevel=5.5

plot file u 1:(($5 )*Yscale-FermiLevel) w l notitle, \
     file u 1:(($6 )*Yscale-FermiLevel) w l notitle, \
     file u 1:(($7 )*Yscale-FermiLevel) w l notitle, \
     file u 1:(($8 )*Yscale-FermiLevel) w l notitle, \
     file u 1:(($9 )*Yscale-FermiLevel) w l notitle, \
     file u 1:(($10)*Yscale-FermiLevel) w l notitle, \
     file u 1:(($11)*Yscale-FermiLevel) w l notitle, \
     file u 1:(($12)*Yscale-FermiLevel) w l notitle, \
     file u 1:(($13)*Yscale-FermiLevel) w l notitle, \
     file u 1:(($14)*Yscale-FermiLevel) w l notitle
