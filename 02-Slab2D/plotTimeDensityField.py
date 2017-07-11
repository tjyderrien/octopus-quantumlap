#!/usr/bin/env python
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2017 T. J.-Y. Derrien
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
from numpy import genfromtxt, loadtxt, chararray
import matplotlib as mp
import matplotlib.pyplot as plt

from matplotlib import rc
# from pylab import *
from scipy.constants import c, epsilon_0, mu_0, pi, e, m_e, h
from matplotlib.legend_handler import HandlerLine2D
import sys

from libAtomicUnits import *
from libPlottingSettings import *

# import calculateVolume as Sample

# Load the n_el and laser files
geometry = "../Si.xyz"
Ne=loadtxt("n_ex", dtype='str', skiprows=9)
laser=loadtxt("laser", dtype='str', skiprows=6)

# get the vector potential from laser file
#TODO: attribute the column number automatically using the header
time = np.asfarray(laser[:,1]) #column 2
A1x = np.asfarray(laser[:,2]) #column number 3
A2x = np.asfarray(laser[:,5]) #column number 6 
timeNe = np.asfarray(Ne[:,1])
density = np.asfarray(Ne[:,2]) #number of particles, not density! 
volume = 1. #Sample.VolumeXYZ(geometry) #SI units! #TODO: get the volume in the octopus log - NOT VALID FOR BULK !!!!
densitySI = density / BohrRadius**3 / volume #Now its a density... TODO: VALIDATE VOLUME CALC.
densityAU = density
# Converting vector potential to field
EtotX = -(np.diff(A1x)/np.diff(time) + np.diff(A2x)/np.diff(time)) / c #E=d_t (A1+A2) / c

# # Convert time from AU to SI
UnitTime = 1E15 #we show in fs
TimeAUtoSI = Time_AU_to_SI(1.) #1.054571800e-34 / 1.6021766208e-19 / 27.21138505e0
time = UnitTime * time * TimeAUtoSI
timeNe = UnitTime * timeNe * TimeAUtoSI
tmax=40e-15*UnitTime

# Plot n_ex and electric field on the same plot
plt.figure()
plt.subplot(211)
#plt.xlabel('Time (fs)')
plt.ylabel('E field (a.u.)')
# plt.xlim((0, time.max()))
plt.xlim((0,tmax))
plt.plot(time[1:], c*EtotX, 'k-', label='c x E(t)')
plt.plot(time, (A1x+A2x), 'b-', label='A(t)')
plt.legend(loc=2)

plt.subplot(212)
plt.xlim((0,time.max()))
plt.plot(timeNe, densityAU, 'k-')  #TODO: change to SI units, and maybe divide by the cell volume
plt.xlabel('Time (fs)')
plt.ylabel('Density (a.u.)')
# plt.legend(loc=2)
plt.savefig('DensityAndField.eps')
plt.savefig('DensityAndField.png')
#plt.show()
