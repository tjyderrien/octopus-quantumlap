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
Ne=loadtxt("../td.general/n_ex", dtype='str', skiprows=9)
#laser=loadtxt("../td.general/laser", dtype='str', skiprows=6) #should be replaced by genfromtxt method, as some errors are present in laser file
laser=np.genfromtxt("../td.general/laser", dtype='str', skip_header=9, loose=False, invalid_raise=True, filling_values=0.)
PostTreatedDensity = loadtxt("postCalculationNexcited.log", dtype='str', skiprows=1)


# get the vector potential from laser file
#TODO: attribute the column number automatically using the header
### LASER PULSE FILE
time = np.asfarray(laser[:,1]) #column 2
A1x = np.asfarray(laser[:,2]) #column number 3
A2x = np.asfarray(laser[:,5]) #column number 6 

timestepTable = np.diff(time) #size N-1
TimeStep = timestepTable[0] #they are all equal, just select the first
A_squared = A1x**2
A_integrated = TimeStep * np.cumsum(A1x) #int_{-\inf}^{t} A(t') dt'
Efield = - (np.diff(A1x)/np.diff(time)) #correct value is *1/c
print np.max(abs(Efield))
 
### N_EX FILE
timeNe = np.asfarray(Ne[:,1])
density = np.asfarray(Ne[:,2]) #number of particles, not density! 

### POST-PROCESSED FILE 
timePost = np.asfarray(PostTreatedDensity[:,2]) #3: time SI
NePost = np.asfarray(PostTreatedDensity[:,3]) #4: N_e
NePostSquaredInside = np.asfarray(PostTreatedDensity[:,6]) #7: N_e squared inside
NePostSquaredOutside = NePost**2 #N_e squared outside the integral

### OTHERS

CellVolume_AU = 0.07820 #Bohr**3
Nproj_AU = density * CellVolume_AU
Npost_AU = NePost * CellVolume_AU
Npost_2in_AU = NePostSquaredInside * CellVolume_AU
Npost_2out_AU= NePostSquaredOutside * CellVolume_AU

# # Convert time from AU to SI
UnitTime = 1E15 #we show in fs
TimeAUtoSI = 1.054571800e-34 / 1.6021766208e-19 / 27.21138505e0
time = UnitTime * time * TimeAUtoSI
timeNe = UnitTime * timeNe * TimeAUtoSI

xmin=0.; xmax=5;
ymax=1E-6; ymax2=1E0

# Plot n_ex and electric field on the same plot
plt.figure()
plt.subplot(211)
#plt.xlabel('Time (fs)')
plt.ylabel('A field (a.u.)')
plt.xlim((0, time.max()))
# plt.xlim((xmin,xmax))
# plt.ylim((-ymax2,ymax2))
# plt.plot(time[1:], c*EtotX, 'k-', label='c x E(t)')
# plt.plot(time, (A1x), 'b-', label='A(t)')
# plt.plot(time, A_integrated, 'r-', label="$\int_{-\infty}^{t} A(t') dt'$")
# plt.plot(time, A_squared, 'g-', label="A^2 (t)")
plt.plot(time[1:], Efield, 'k-', label="|E(t)|") 
plt.legend(loc=2)

plt.subplot(212)
plt.xlim((0,time.max()))
# plt.xlim((xmin,xmax))
# plt.ylim((0,ymax))
plt.plot(timeNe, Nproj_AU, 'k-', label=r'$N_e^{proj}$')  #TODO: change to SI units, and maybe divide by the cell volume
plt.plot(timeNe, Npost_AU, 'b-', label=r'$\int_{\infty}^{t} \left( N_{td}(t\'{}) - N_{gs} \right) dt\'{}$')
plt.plot(timeNe, 100.*Npost_2in_AU, 'r-', label=r'$\int_{\infty}^{t} \left( N_{td}(t\'{}) - N_{gs} \right)^2 dt\'{}$')
plt.plot(timeNe, 1.6/23. * Npost_2out_AU, 'g-', label=r'$\left( \int_{\infty}^{t} \left( N_{td}(t\'{}) - N_{gs} \right) dt\'{} \right)^2 $')
plt.xlabel('Time (fs)')
plt.ylabel('Density (a.u.)')
plt.legend(loc=2)
plt.savefig('KeldyshAnalysis.eps')
plt.savefig('KeldyshAnalysis.png')
plt.show()
