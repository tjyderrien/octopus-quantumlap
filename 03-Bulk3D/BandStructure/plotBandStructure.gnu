#!gnuplot

# set terminal x11

file="bandstructure"

#Column meanings. 
# 1 2  3  4  5     6     7     ... 8 9 
# ? kx ky kz band1 band2 band3 ....

plot file u 1:5 w lp notitle, \
     file u 1:6 w lp notitle, \
     file u 1:7 w lp notitle, \
     file u 1:8 w lp notitle, \
     file u 1:9 w lp notitle, \
     file u 1:10 w lp notitle, \
     file u 1:11 w lp notitle, \
     file u 1:12 w lp notitle, \
     file u 1:13 w lp notitle, \
     file u 1:14 w lp notitle
