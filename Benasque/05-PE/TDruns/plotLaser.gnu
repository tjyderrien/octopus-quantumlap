#!gnuplot

# set terminal dumb
set isosamples 200
set samples 200

A(t)=A0*cos(pi/2. * (t - t0 - 2.*tau0)/tau0)*cos(omega*t)
A0=0.5
omega=2.214088
tau=2.*pi/omega
t0=10.*tau
tau0=t0

#set xrange [-100:100]
plot A(x) w l t 'gnuplot', \
"td.general/laser" u 2:3 w l lc 3 t 'octopus'
