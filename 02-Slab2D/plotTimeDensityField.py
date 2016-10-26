#!/usr/bin/env python
#-*- coding: utf-8 -*-

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

import calculateVolume as Sample

# Settings for matplotlib
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':'16'})
## for Palatino and other serif fonts use:
#rc('font', **{'family':'serif', 'serif':['Palatino'], 'size':'18'})
#rc('text', usetex=True)
mp.rcParams['legend.numpoints'] = 1

# Some constants
BohrRadius=5.2917721067e-11 #m

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
volume = Sample.VolumeXYZ(geometry)
densitySI = density / BohrRadius**3 / volume #Now its a density... TODO: VALIDATE VOLUME CALC.
# Convert vector potential to field
EtotX = -(np.diff(A1x)/np.diff(time) + np.diff(A2x)/np.diff(time)) / c #E=d_t (A1+A2) / c

# # Convert time from AU to SI
UnitTime = 1E15 #we show in fs
TimeAUtoSI = 1.054571800e-34 / 1.6021766208e-19 / 27.21138505e0
time = UnitTime * time * TimeAUtoSI
timeNe = UnitTime * timeNe * TimeAUtoSI

# Plot n_ex and electric field on the same plot
plt.figure()
plt.subplot(211)
#plt.xlabel('Time (fs)')
plt.ylabel('E field (a.u.)')
plt.xlim((0, time.max()))
plt.plot(time[1:], c*EtotX, 'k-', label='c x E(t)')
plt.plot(time, (A1x+A2x), 'b-', label='A(t)')
plt.legend(loc=2)
plt.subplot(212)
plt.xlim((0,time.max()))
plt.plot(timeNe, densitySI, 'k-')  #TODO: change to SI units, and maybe divide by the cell volume
plt.xlabel('Time (fs)')
plt.ylabel('Density (a.u.)')
# plt.legend(loc=2)
plt.savefig('DensityAndField.eps')
plt.savefig('DensityAndField.png')
#plt.show()
