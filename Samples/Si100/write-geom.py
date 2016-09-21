#!/usr/bin/python

import sys
from math import sqrt

A = 5.431 #5.431: exp @ 300 K, 5.430: exp @ 0K, 5.416: LDA-value ?
L_vac = 50.
D = A/sqrt(2.)
if len(sys.argv) > 1:
    n = int(sys.argv[1]);
    m = int(sys.argv[2]);
#    L = int(sys.argv[3]);
else:
    n = 1
    m = 1 
#    L = 30.
# lattice vectors
print "lattice_vector  %.16f  %.16f  %.16f" % (n*A+L_vac, 0., 0.) #x
print "lattice_vector  %.16f  %.16f  %.16f" % (0, m*D, 0.)  #y
print "lattice_vector  %.16f  %.16f  %.16f" % (0, 0., m*D)  #z

C = A / (2.*sqrt(2.))
Bohr = 0.52918 #Angstrums

# Vectors for adding the H atoms
# U1x = (4.3884-2.1942)*Bohr; iU1z=-2.1942*Bohr #/!\ Bohrs
vecA = 1.81394*Bohr; vecB = 4.3884*Bohr; vecC = 2.1942*Bohr; vecD = -1.73952*Bohr
#FrontSide
U1x = -vecD; U1y = -vecC;        U1z = 0.;
U2x = -vecD; U2y =  vecB - vecC; U2z = 0.;
#Backside
U3x =  vecD; U3y = -vecC;        U3z = 0.;
U4x =  vecD; U4y =  vecB - vecC; U4z = 0.;

# backside layer
# for j in range(m):
#    for k in range(m):
#        print "atom %.16f  %.16f  %.16f  H" % (0*A+U3x, D*(k)+U3y, D*(j)+U3z)
#        print "atom %.16f  %.16f  %.16f  H" % (0*A+U4x, D*(k)+U4y, D*(j)+U4z)

# generate the lattice positions
for i in range(n): #careful: X is depth here. Backside: X=0. Surface: X=A*n.
    for j in range(m): #Y: we keep same distance in directions // to surface. ExtremeLeft = -A*m/2. ExtremeRight= A*m/2.
	for k in range(m): #Z: see Y. 
		print "atom %.16f  %.16f  %.16f  Si" % (i*A, D*(k), D*(j))
		print "atom %.16f  %.16f  %.16f  Si" % ((i+0.25)*A, D*(k), C+D*(j))
		print "atom %.16f  %.16f  %.16f  Si" % ((i+0.5)*A, C+D*(k),  C+D*(j))
		print "atom %.16f  %.16f  %.16f  Si" % ((i+0.75)*A, C+D*(k),  D*(j))

#front side
#for j in range(m):
#    for k in range(m):
#        print "atom %.16f  %.16f  %.16f  H" % ((i+0.75)*A+U1x, C+D*(k)+U1y, D*(j)+U1z)
#        print "atom %.16f  %.16f  %.16f  H" % ((i+0.75)*A+U2x, C+D*(k)+U2y, D*(j)+U2z)
