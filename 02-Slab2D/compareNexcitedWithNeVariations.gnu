#!gnuplot

reset

set xlabel 'Time (fs)' 
set ylabel 'Density (at.u.)'

set output 'Nexcited.eps'
set terminal postscript eps enhanced color font 'Helvetica, 24'

ProjectionFile="../td.general/n_ex"
PostTreatFile="postCalculationNexcited.log"

AUtoSI_time="2.44e-17"

plot PostTreatFile u 3:4 w l lw 3 lc 1 t 'Ne(t) post-treatment', \
ProjectionFile u ($2*AUtoSI_time):3 w l t 'N_{ex}(t)=< {/Symbol f}| {Symbol f}_{GS} >' 
