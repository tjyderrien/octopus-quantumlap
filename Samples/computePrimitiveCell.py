#!/usr/bin/env python2
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2020 T. J.-Y. Derrien
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>


## @package libPrimitiveCells
# Defines primitive cells 
# These routines already exist in softwares like Avogadro, or library for group theory, or ASE, but somehow I start to doubt them. 
# See the 7 cases of the Kittel book. 

import numpy as np

angstrom=1
deg=np.pi/180.


#a=3.10321 #angstrom
#alpha = 109.47122 #deg

## Cubic symmetric crystalline structure
# a=b=c and alpha=beta=gamma. 
# We have 3 vectors, a, b, c. 
# a=<a, 0, 0>. 
# b=<b1, b2, 0>
# c=<c1, c2, c3>
def cubic(a, alpha):
    a=a*angstrom; alpha=alpha*deg
    b=a; c=a;
    beta=alpha; gamma=alpha;
    a1=a; a2=0; a3=0; b3 = 0
    
    # Branch selection
    s1 = 1 #or -1
    
    b1 = a*np.cos(gamma)
    b2 = s1*a*np.sqrt(1-np.cos(gamma)**2)
    
    c1 = a*np.cos(beta)
    c2 = (a**2*np.cos(alpha)-(a*np.cos(gamma))*(a*np.cos(beta)))/b2
    c3 = np.sqrt(a**2 - c1**2-c2**2)
    
    print "Not normalized"
    print str(a1)+"|"+str(a2)+"|"+str(a3)
    print str(b1)+"|"+str(b2)+"|"+str(b3)
    print str(c1)+"|"+str(c2)+"|"+str(c3)
    print ""
    print "Normalized"
    print str(a1/a)+"|"+str(a2/a)+"|"+str(a3/a)
    print str(b1/a)+"|"+str(b2/a)+"|"+str(b3/a)
    print str(c1/a)+"|"+str(c2/a)+"|"+str(c3/a)
  

## Cubic symmetric crystalline structure
# a=b !=c, alpha=beta != gamma
# We have 3 vectors, a, b, c. 
# a=<a, 0, 0>. 
# b=<b1, b2, 0>
# c=<c1, c2, c3>
def hexagonal(a, c, alpha, gamma):
    a=a*angstrom; alpha=alpha*deg
    b=a; 
    c=c*angstrom; gamma=gamma*deg
    beta=alpha; 
    a1=a; a2=0; a3=0; b3 = 0
    
    # Branch selection for b2
    s1 = 1 #or -1
    
    b1 = a*np.cos(gamma)
    c1 = c*np.cos(alpha)
    
    b2 = s1*np.sqrt(a**2-b1**2)
    
    c2 = (a*c*np.cos(alpha)-b1*c1)/b2
    c3 = s2*np.sqrt(c**2 - c1**2-c2**2)
    
    print "Not normalized"
    print str(a1)+"|"+str(a2)+"|"+str(a3)
    print str(b1)+"|"+str(b2)+"|"+str(b3)
    print str(c1)+"|"+str(c2)+"|"+str(c3)
    print ""
    print "Normalized"
    print str(a1/a)+"|"+str(a2/a)+"|"+str(a3/a)
    print str(b1/a)+"|"+str(b2/a)+"|"+str(b3/a)
    print str(c1/c)+"|"+str(c2/c)+"|"+str(c3/c)
