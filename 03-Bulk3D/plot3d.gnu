#!gnuplot

reset
summary="Database.csv"

## PLOT 1: N_final(Efield, CEP1)
# Filtering single pulses
splot "< awk '{ if($6<1.) print }' Database.csv" u 3:4:13 w p
