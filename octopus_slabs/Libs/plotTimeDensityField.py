#!/usr/bin/env python
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2023 T. J.-Y. Derrien
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


# IMPORT LIBRARIES
import numpy as np
from numpy import loadtxt
# from pylab import *
from scipy.constants import c
import matplotlib.pyplot as plt

from octopus_slabs.Libs.libAtomicUnits import *

# import calculateVolume as Sample

# Load the n_el and laser files
geometry = "../Si.xyz"
# Ne=loadtxt("n_ex", dtype='str', skiprows=9)
laser=loadtxt("laser.tmp", dtype='str', skiprows=6)

# get the vector potential from laser file
#TODO: attribute the column number automatically using the header
time = np.asfarray(laser[:,1]) #column 2
A1x = np.asfarray(laser[:,2]) #column number 3
A2x = np.asfarray(laser[:,5]) #column number 6 
# timeNe = np.asfarray(Ne[:,1])
# density = np.asfarray(Ne[:,2]) #number of particles, not density! 
volume = 1. #Sample.VolumeXYZ(geometry) #SI units! #TODO: get the volume in the octopus log - NOT VALID FOR BULK !!!!
# densitySI = density / BohrRadius**3 / volume #Now its a density... TODO: VALIDATE VOLUME CALC.
# densityAU = density
# Converting vector potential to field

E1x = -np.gradient(A1x)/np.gradient(time) / c
E2x = -np.gradient(A2x)/np.gradient(time) / c

EtotX = E1x+E2x


# # Convert time from AU to SI
UnitTime = 1E15 #we show in fs
TimeAUtoSI = Time_AU_to_SI(1.) #1.054571800e-34 / 1.6021766208e-19 / 27.21138505e0
time = UnitTime * time * TimeAUtoSI
# timeNe = UnitTime * timeNe * TimeAUtoSI

tmax=100e-15*UnitTime

# Plot n_ex and electric field on the same plot
plt.figure()
ax1 = plt.subplot(111)
ax11 = ax1.twinx()

ax1.set_xlabel('Time (fs)')
ax1.set_ylabel('E field (a.u.)')
ax11.set_ylabel('E field (a.u.)')
# ax1.set_ylabel('A field (a.u.)')
# plt.xlim((0, time.max()))
ax1.set_xlim((0,tmax))
ax11.set_xlim((0,tmax))
# plt.plot(time[1:], c*EtotX, 'k-', label='c x E(t)')
ax1.plot(time, E1x, 'b.', label=r'$E_1(t)$')
ax11.plot(time, E2x, 'r--', label=r'$E_2(t)$')
# ax11.plot(time, A2x, 'r-', label=r'$A_2(t)$')
# plt.legend(loc=2)


"""plt.subplot(212)
plt.xlim((0,time.max()))
plt.plot(timeNe, densityAU, 'k-')  #TODO: change to SI units, and maybe divide by the cell volume
plt.xlabel('Time (fs)')
plt.ylabel('Density (a.u.)')
"""
# plt.legend(loc=2)

plt.savefig('Field.eps')
plt.savefig('Field.png')
# plt.show()
