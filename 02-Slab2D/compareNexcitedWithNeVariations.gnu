#!gnuplot

reset

set xlabel 'Time (fs)' 
set ylabel 'Density [number of e-]'

unitx=1E15
outputname='Nexcited.eps'

set output outputname
set terminal postscript eps enhanced color font 'Helvetica, 24'
# set terminal dumb

# set xrange [:5e-15*unitx]
set key out bottom center
unset log y

ProjectionFile="../td.general/n_ex"
PostTreatFile="postCalculationNexcited.log"
CellVolumeAU=system(" grep 'Cell volume' $( ls -1htr ../*.log | tail -n1 ) | awk '{ print $4 }' ") #Bohr^3
AUtoSI_time="2.44e-17"
Bohr=5.2917721067e-11 #m
CellVolume=CellVolumeAU*(Bohr**3)

# print ""
print "Volume=", CellVolumeAU, "[Bohr**3] = ", CellVolume, "[m^-3]."

plot PostTreatFile u ($3*unitx):($4) w l lw 3 lc 1 t 'N_{e}(t) post-treatment', \
PostTreatFile u ($3*unitx):($5) w l lw 3 lc 3 t 'N_{h}(t) post-treatment', \
PostTreatFile u ($3*unitx):(($4+$5)) w l lw 3 lc 2 t '(N_e+N_{h})(t) post-treatment', \
ProjectionFile u ($2*AUtoSI_time*unitx):3 w l lc 4 lw 2 t 'N_{ex}(t)=< {/Symbol f} | {/Symbol f}_{GS} >' 

print "EPS file was written."

set output
set terminal dumb
replot

print "** Info: ".outputname." EPS file was also written."
