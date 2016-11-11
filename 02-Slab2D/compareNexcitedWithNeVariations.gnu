#!gnuplot

reset
pwd
set xlabel 'Time (fs)' 
set ylabel 'Number of e-'
set y2label 'Total number of e-'

unitx=1E15
outputname='Nexcited.eps'

set output outputname
set terminal postscript eps enhanced color font 'Helvetica, 24'
# set terminal dumb


set xrange [:40e-15*unitx]
set key out top center

unset log y

set ytics nomirror
set y2tics mirror


ProjectionFile="../td.general/n_ex"
PostTreatFile="postCalculationNexcited.log"
UnitCellVolumeAU=system(" grep 'volume/point' $( ls -1htr ../*.log | tail -n1 ) | awk '{ print $11 }' ") #Bohr^3
AUtoSI_time="2.418884326505e-17" #https://en.wikipedia.org/wiki/Atomic_units
Bohr=5.2917721067e-11 #m
CellVolume=UnitCellVolumeAU*(Bohr**3) #Volume of one element

# print ""
print "Element volume=", UnitCellVolumeAU, "[Bohr**3] = ", CellVolume, "[m^-3]."
print "Spin: |< psi_k (t) | psi_k (GS) >|^2 gives density for ONE spin only."

plot PostTreatFile u ($3*unitx):(8.17/1.34*$4*UnitCellVolumeAU) w l lw 3 lc 1 t 'N_{e}(t) post-treatment', \
ProjectionFile u ($2*AUtoSI_time*unitx):($3) w l lc 2 lw 2 t 'N_{ex}(t)=< {/Symbol f} | {/Symbol f}_{GS} >', \
PostTreatFile u ($3*unitx):($7*UnitCellVolumeAU) w l lw 3 lc 4 t 'int[(TD-GS)**2]', \
PostTreatFile u ($3*unitx):(8.16/1.74*($4*UnitCellVolumeAU)**2) w l lw 3 lc 5 t 'int[(TD-GS)]**2'

# PostTreatFile u ($3*unitx):(($4+$5)*UnitCellVolumeAU) w l lw 3 lc 2 t '(N_e+N_{h})(t) post-treatment', \
# PostTreatFile u ($3*unitx):(-$5*UnitCellVolumeAU) w l lw 3 lc 3 t 'N_{h}(t) post-treatment', \
# PostTreatFile u ($3*unitx):($6*UnitCellVolumeAU) w l lw 3 lc 3 axis x1y2 t 'N_{tot}', \

print "** Info: ".outputname." EPS file was also written."

set output
set terminal dumb
replot

# print "** Info: ".outputname." EPS file was also written."
