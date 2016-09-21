#!/usr/bin/python

import sys
from math import sqrt

A = 5.416
L_vac = 30.
D = A/sqrt(2.)
if len(sys.argv) > 1:
    n = int(sys.argv[1]) #number of cells normal to surface
    p = int(sys.argv[2]) #number of cells parallel to surface
else:
    n = 1
    p = 1

# lattice vectors
print "lattice_vector  %.16f  %.16f  %.16f" % (A/sqrt(2.), 0., 0.) #x
print "lattice_vector  %.16f  %.16f  %.16f" % (0., A/sqrt(2.), 0.) #y
print "lattice_vector  %.16f  %.16f  %.16f" % (0., 0., n*A+L_vac) #z

for i in range(n):
	for j in range(p): 
		for k in range(p): 
			C = A / (2.*sqrt(2.))
			print "atom %.16f  %.16f  %.16f  Si" % (D*j, D*k, i*A)
			print "atom %.16f  %.16f  %.16f  Si" % (D*j, C+D*k,  (i+0.25)*A)
			print "atom %.16f  %.16f  %.16f  Si" % (C+D*j,  C+D*k,  (i+0.5)*A)
			print "atom %.16f  %.16f  %.16f  Si" % (C+D*j,  D*k, (i+0.75)*A)
