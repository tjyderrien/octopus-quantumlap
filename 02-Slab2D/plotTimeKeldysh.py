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

TotalCurrent=loadtxt("../td.general/total_current", dtype='str', skiprows=4)


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

### TOTAL_CURRENT file
ShrinkBy=3
timeCurrent=np.asfarray(TotalCurrent[:,1])
timeCurrent=timeCurrent[:timeCurrent.size-ShrinkBy]
CurrentX=np.asfarray(TotalCurrent[:,2])
CurrentX=CurrentX[:CurrentX.size-ShrinkBy]

### POST-PROCESSED FILE 
timePost = np.asfarray(PostTreatedDensity[:,2]) #3: time SI
NePost = np.asfarray(PostTreatedDensity[:,3]) #4: N_e
NePostSquaredInside = np.asfarray(PostTreatedDensity[:,6]) #7: N_e squared inside
NePostSquaredOutside = NePost**2 #N_e squared outside the integral

### OTHERS

CellVolume_AU = Sample.VolumeLog("../") #0.06385 #0.07820 #Bohr**3
print "Cell volume = "+str(CellVolume_AU)+" Bohr**3"
Nproj_AU = density #* CellVolume_AU
Npost_AU = NePost * CellVolume_AU
Npost_2in_AU = NePostSquaredInside * CellVolume_AU
Npost_2out_AU= NePostSquaredOutside * CellVolume_AU
NnewDerivDt=np.multiply(CurrentX,A1x) #still to multiply by dt and sum.
Nnew = np.cumsum(NnewDerivDt)

# # Convert time from AU to SI
UnitTime = 1E15 #we show in fs
TimeAUtoSI = 1.054571800e-34 / 1.6021766208e-19 / 27.21138505e0
time = UnitTime * time * TimeAUtoSI
timeNe = UnitTime * timeNe * TimeAUtoSI
timeCurrentdt = timeCurrent #* TimeAUtoSI
timeCurrent=UnitTime * timeCurrent * TimeAUtoSI
dtCurrent = timeCurrentdt[1]-timeCurrentdt[0] 
print "TimeStep = "+str(dtCurrent)

xmin=0.; xmax=40.; unitx=1E15
ymax=1E-6; ymax2=1E0

# Plot n_ex and electric field on the same plot
plt.figure()
plt.subplot(211)
#plt.xlabel('Time (fs)')
plt.ylabel('A field (a.u.)')
# plt.xlim((0, time.max()))
plt.xlim((xmin,xmax))
# plt.ylim((-ymax2,ymax2))
# plt.plot(time[1:], c*EtotX, 'k-', label='c x E(t)')
plt.plot(time, (A1x), 'b-', label='A(t)')
# plt.plot(time, A_integrated, 'r-', label="$\int_{-\infty}^{t} A(t') dt'$")
# plt.plot(time, A_squared, 'g-', label="A^2 (t)")
plt.plot(time[1:], Efield, 'k-', label="|E(t)|")
wavelength_SI = 800e-9 #TODO: make automatic!
wavelength_AU = wavelength_SI/BohrRadius
m=np.arange(0,100) #[-5.,-4.,-3.,-2.,-1.,0.,1.,2.,3.,4.,5.]
xposition = np.multiply(m,unitx*wavelength_SI/(2.*c))-3.*pi/2.
for xc in xposition:
  plt.axvline(x=xc, color='k', linestyle='--')
plt.legend(loc=2)

plt.subplot(212)
# plt.xlim((0,unitx*timePost.max()))
plt.xlim((xmin,xmax))
# plt.ylim((0,ymax))
plt.plot(timeNe, Nproj_AU, 'k-', label=r'$N_e^{proj}(t)$')  #TODO: change to SI units, and maybe divide by the cell volume
plt.plot(unitx*timePost, Npost_AU, 'b-', label=r'$\int_{\Omega} \left( n_{td}(r,t) - n_{gs}(r) \right)^{+} d^{3}r$')
for xc in xposition:
  plt.axvline(x=xc, color='k', linestyle='--')
# plt.plot(unitx*timePost, 10.*Npost_2in_AU, 'r-', label=r'$\int_{\Omega} \left( n_{td}(r,t)) - n_{gs}(r) \right)^2 d^{r}$')
# plt.plot(unitx*timePost, 1./80. * Npost_2out_AU, 'g-', label=r'$\left[ \int_{\Omega} ( N_{td}(r,t) - N_{gs}(r) ) d^{3}r \right]^{2}$')
print timeCurrent.shape, Nnew.shape
# plt.plot(timeCurrent, Nnew/10.26**3, 'g-', label=r'$\int_{-\infty}^{t} J(t^{\prime}) A(t^\prime) dt^{\prime} $') #TODO: still a problem in units
plt.xlabel('Time (fs)')
plt.ylabel('Density (a.u.)')
# plt.legend(loc=2)
plt.legend(bbox_to_anchor=(0, -0.6, 1, 0), loc='right', ncol=2, mode="expand", borderaxespad=0)
plt.savefig('KeldyshAnalysis.eps', bbox_inches='tight')
plt.savefig('KeldyshAnalysis.png', bbox_inches='tight')
# plt.show()
