#!gnuplot

reset
summary="Database.csv"

plot1="< awk '{ if($6<1.) print }'".summary

## PLOT 1: N_final(Efield, CEP1)
splot "< awk '{ if($6<1.) print }' Database.csv" u 3:4:13 w p
