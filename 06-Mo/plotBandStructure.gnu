#!gnuplot

# set terminal x11

set terminal postscript eps enhanced font 'Helvetica, 24'
set output 'BandStructure.eps'


file="bandstructure"

HartreeToEv = 27.2114 #1 Hartree = 27.2114 eV

Yscale = HartreeToEv

set ylabel 'Energy (eV)'
set yrange [-10:10]
unset xtics

#Column meanings. 
# 1 2  3  4  5     6     7     ... 8 9 
# ? kx ky kz band1 band2 band3 ....


FermiLevel=0.451415*HartreeToEv

plot file u 1:(($5 )*Yscale-FermiLevel) w l notitle, \
     file u 1:(($6 )*Yscale-FermiLevel) w l notitle, \
     file u 1:(($7 )*Yscale-FermiLevel) w l notitle, \
     file u 1:(($8 )*Yscale-FermiLevel) w l notitle, \
     file u 1:(($9 )*Yscale-FermiLevel) w l notitle, \
     file u 1:(($10)*Yscale-FermiLevel) w l notitle, \
     file u 1:(($11)*Yscale-FermiLevel) w l notitle, \
     file u 1:(($12)*Yscale-FermiLevel) w l notitle, \
     file u 1:(($13)*Yscale-FermiLevel) w l notitle, \
     file u 1:(($14)*Yscale-FermiLevel) w l notitle, \
     file u 1:(($15)*Yscale-FermiLevel) w l notitle, \
     file u 1:(($16)*Yscale-FermiLevel) w l notitle, \
     file u 1:(($17)*Yscale-FermiLevel) w l notitle, \
     file u 1:(($18)*Yscale-FermiLevel) w l notitle, \
     file u 1:(($19)*Yscale-FermiLevel) w l notitle, \
     file u 1:(($20)*Yscale-FermiLevel) w l notitle, \
     file u 1:(($21)*Yscale-FermiLevel) w l notitle, \
     file u 1:(($22)*Yscale-FermiLevel) w l notitle, \
     file u 1:(($23)*Yscale-FermiLevel) w l notitle, \
     file u 1:(($24)*Yscale-FermiLevel) w l notitle

#   0    | 0     | 0     #Gamma
#   0.5  | -0.5  | 0.5   #H
#   0    | 0     | 0.5   #N
#   0    | 0     | 0.0   #Gamma
#   0.25 | 0.25  | 0.25  #P
#   0.5  | -0.5  | 0.5   #H(|P)
#   0    | 0     | 0.5   #N







