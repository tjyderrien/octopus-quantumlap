#!/usr/bin/python

import sys
from math import sqrt

A = 5.416
L_vac = 30.

if len(sys.argv) > 1:
    n = int(sys.argv[1]);
    m = int(sys.argv[2]);
else:
    n = 1
    m = 1 
# lattice vectors
print("lattice_vector  %.16f  %.16f  %.16f" % (n*A+L_vac, 0., 0.)) #x
print("lattice_vector  %.16f  %.16f  %.16f" % (0, A/sqrt(2), 0.))  #y
print("lattice_vector  %.16f  %.16f  %.16f" % (0, 0., A/sqrt(2)))  #z

# generate the lattice positions
for i in range(n): #careful: X is depth here. 
    for j in range(m): #we keep same distance in directions // to surface
	C = A / (2*sqrt(2.))
	print("atom %.16f  %.16f  %.16f  Si" % (i*A, 0., 0.))
	print("atom %.16f  %.16f  %.16f  Si" % ((i+0.25)*A, 0., C))
	print("atom %.16f  %.16f  %.16f  Si" % ((i+0.5)*A, C,  C))
	print("atom %.16f  %.16f  %.16f  Si" % ((i+0.75)*A, C,  0.))
