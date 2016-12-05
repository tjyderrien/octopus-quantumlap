#!/usr/bin/env python
#-*- coding: utf-8 -*-

# IMPORT LIBRARIES
import numpy as np
from numpy import genfromtxt, loadtxt, chararray
# from pylab import *
import sys

import calculateVolume as Sample
import libKeldysh as Keldysh
from libAtomicUnits import *
from libPlottingSettings import *
from libString import *

# TODO: choose AU or SI from beginning, and centralize their calculation
# TODO: separate plotting and density post-calculation into libs

UnitSystem="SI" #Defines the system unit for the observables. Time will remain in SI. 
QuantityOfInterest = "density" #Would you prefer to work densities, or particle number? 

# Load the n_el and laser files (atomic units)
geometry = "../Si.xyz"
Ne=loadtxt("../td.general/n_ex", dtype='str', skiprows=9) #total number of excited electrons, not a density!
#laser=loadtxt("../td.general/laser", dtype='str', skiprows=6) #should be replaced by genfromtxt method, as some errors are present in laser file
laser=np.genfromtxt("../td.general/laser", dtype='str', skip_header=9, loose=False, invalid_raise=True, filling_values=0.)
PostTreatedDensity = loadtxt("postCalculationNexcited.log", dtype='str', skiprows=1) #Unit: number of electrons, which should be multiplied by cell volume (AU). 
TotalCurrent=loadtxt("../td.general/total_current", dtype='str', skiprows=4)

# get the vector potential from laser file
#TODO: attribute the column number automatically using the header
### LASER PULSE FILE
time = np.asfarray(laser[:,1]) #column 2
A1x = np.asfarray(laser[:,2]) #column number 3
A2x = np.asfarray(laser[:,5]) #column number 6 

c_AU = 137.

timestepTable = np.diff(time) #size N-1
TimeStep = timestepTable[0] #they are all equal, just select the first
A_squared = A1x**2
A_integrated = TimeStep * np.cumsum(A1x) #int_{-\inf}^{t} A(t') dt'
Efield = - (np.diff(A1x)/np.diff(time)) #correct value is *1/c 
Efield_amplitude_SI=np.max(abs(Field_AU_to_SI(Efield / c_AU)))
print "** Info: Electric field amplitude = "+str(Efield_amplitude_SI / 1E9)+" V/nm."
 
### N_Keldysh and N_Gruzdev
# Capturing data from input file
wavelength_inp = LookUpStringInFile("../inp", "lambdaA =")
wavelength_SI = float(wavelength_inp[10:17])
PulseDuration_inp = LookUpStringInFile("../inp", "tauAsi = ")
PulseDuration_SI = float(PulseDuration_inp[8:15])

# Explicit data
# PulseDuration_SI = 10E-15
Egap_SI = 2.58*e

meff = 0.18
# meff = 0.43
# wavelength_SI = 800e-9
ValenceBandMax = 5E28 #TODO: HOW TO INTRODUCE THIS in Keldysh?
dtPrecision = 1E-17
order=50
CentralTime = 2.*PulseDuration_SI #20E-15
ShowPlot=False

PeakFluence_SI = Keldysh.FieldToIntensity(Efield_amplitude_SI) * PulseDuration_SI

timeKeldysh, N_Keldysh_SI, N_Gruzdev_SI = Keldysh.plotPulseToDensity(Egap_SI, meff, wavelength_SI, PulseDuration_SI, PeakFluence_SI, dtPrecision, order, ShowPlot, CentralTime)  

### N_EX octopus FILE
timeNe = np.asfarray(Ne[:,1])
Nproj_AU = np.asfarray(Ne[:,2]) #number of particles, not a density!

### TOTAL_CURRENT file #Just for futur investigations with current
ShrinkBy=3
timeCurrent=np.asfarray(TotalCurrent[:,1])
timeCurrent=timeCurrent[:timeCurrent.size-ShrinkBy]
CurrentX=np.asfarray(TotalCurrent[:,2])
CurrentX=CurrentX[:CurrentX.size-ShrinkBy]

### POST-PROCESSED FILE #Some attemps with charge density
timePost = np.asfarray(PostTreatedDensity[:,2]) #3: time SI
NePost = np.asfarray(PostTreatedDensity[:,3]) #4: N_e number of electrons to multiply by cell volume (AU)
NePostSquaredInside = np.asfarray(PostTreatedDensity[:,6]) #7: N_e squared inside
NePostSquaredOutside = NePost**2 #N_e squared outside the integral

### SingleCell and Supercell Volume calculation (atomic units)
SuperCellVolume_AU = Sample.VolumeOfSupercell_Log("../")
SuperCellVolume_SI = Volume_AU_to_SI(SuperCellVolume_AU)
SampleCellVolume_AU = SuperCellVolume_AU - Sample.VolumeOfVacuum("../") #Volume of the sample (without vacuum region) 
# SampleCellVolume_AU = Sample.VolumeXYZ("../") ##Also wrong apparently...
SampleCellVolume_SI = Volume_AU_to_SI(SampleCellVolume_AU)
CellVolume_AU = Sample.VolumePerPoint_Log("../") #0.06385 #0.07820 #Bohr**3
CellVolume_SI = Volume_AU_to_SI(CellVolume_AU)

print ""
print "Super-cell total volume = "+str(SuperCellVolume_AU)+" Bohr**3."
print "Super-cell total volume = "+str(SuperCellVolume_SI)+" m^{3}." 
print ""
print "Sample volume = "+str(SampleCellVolume_AU)+" Bohr**3."
print "Sample vacuum volume = "+str(SampleCellVolume_SI)+" m^{3}."
print ""
print "Cell volume = "+str(CellVolume_AU)+" Bohr**3."
print "Cell volume = "+str(Volume_AU_to_SI(CellVolume_AU))+" m^{-3}."

### Integration in space

# Calculation using quantum projections < phi_gs | phi(t) >
#Nproj_AU = density #* CellVolume_AU 
# Nproj_AU = Nproj_AU * CellVolume_AU #Should we really multiply by volume ?
density_proj_AU = Nproj_AU / SampleCellVolume_AU #This is now a density (a.u.)
# Calculation of particle number from charge density (octopus.output.density)
Npost_AU = NePost * CellVolume_AU #NOW, it's a number of electrons. :)
Npost_2in_AU = NePostSquaredInside * CellVolume_AU #Same
Npost_2out_AU= NePostSquaredOutside * CellVolume_AU #Same

density_post_AU = Npost_AU / SuperCellVolume_AU

# print CurrentX.size, A1x.size
# Shrink A1x to size of CurrentX

### Attempt using J(r,t) x A(r,t). #Purpose? 
Overlap=np.array([CurrentX.size, A1x.size])
indexmax=Overlap.min()
A1xShrinked = A1x[0:indexmax]
CurrentXShrinked=A1x[0:indexmax]
NnewDerivDt=np.multiply(CurrentXShrinked,A1xShrinked) #still to multiply by dt and sum.
Nnew = np.cumsum(NnewDerivDt)

### Convert time from AU to SI
UnitTime = 1E15 #we show time in fs
TimeAUtoSI = Time_AU_to_SI(1.)
time = UnitTime * time * TimeAUtoSI
timeNe = UnitTime * timeNe * TimeAUtoSI

timeCurrentdt = timeCurrent #* TimeAUtoSI
timeCurrent=UnitTime * timeCurrent * TimeAUtoSI
dtCurrent = timeCurrentdt[1]-timeCurrentdt[0]

### Convert plotted quantities from AU to SI
# Avoid ill-defined units
if(UnitSystem == "SI"):
  print "** Info: SI system chosed."
  observable_proj = Density_AU_to_SI(density_proj_AU)
  observable_post = Density_AU_to_SI(density_post_AU)
elif(UnitSystem == "AU"):
  print "** Info: AU system chosed."
  observable_proj = density_proj_AU
  observable_post = density_post_AU

else: 
  print "** Error: Unit system was ill defined. Set UnitSystem=[SI,AU] in plotTimeKeldysh.py."
  exit()
### Convert density from AU to SI
#Density_AU_to_SI()
# It is actually done on the plot, directly. 


### Convert fields from AU to SI

print "TimeStep = "+str(dtCurrent)

### Defining plotting boundaries
xmin=0.; xmax=40.; unitx=1E15
ymax=1E-6; ymax2=1E0

# Plotting electric field
plt.figure(figsize=(10,10))
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
# wavelength_SI is now automatic
wavelength_AU = length_SI_to_AU(wavelength_SI)

#Let's add a grid to compare the extrema
m=np.arange(0,100) #[-5.,-4.,-3.,-2.,-1.,0.,1.,2.,3.,4.,5.]
xposition = np.multiply(m,unitx*wavelength_SI/(2.*c))-3.*pi/2.
for xc in xposition:
  plt.axvline(x=xc, color='k', linestyle='--')
plt.legend(loc=2)

### Plotting density evolution
plt.subplot(212)
# plt.xlim((0,unitx*timePost.max()))
plt.xlim((xmin,xmax))
# plt.ylim((0,ymax))

if(QuantityOfInterest == "density"):
  if(UnitSystem=="AU"):
    plt.ylabel('Density exc. el. (a.u.)')
  elif(UnitSystem=="SI"):
    plt.ylabel('Density exc. el. (m^{-3})')
  plt.plot(timeNe, observable_proj, 'k-', label=r'TD-DFT: $n_e^{proj}(t)$')
#  plt.plot(unitx*timePost, observable_post, 'k--', label=r'TD-DFT: \frac{1}{V_{tot}} $\int_{\Omega} \left( n_{td}(r,t) - n_{gs}(r) \right)^{+} d^{3}r$')
  plt.plot(timeKeldysh*unitx, N_Keldysh_SI, 'r-', label=r'[Keldysh (1964)]')
  plt.plot(timeKeldysh*unitx, N_Gruzdev_SI, 'b-', label=r'[Gruzdev (2004)]')
elif (QuantityOfInterest == "number"):
  plt.ylabel('Number of exc. electrons')
  plt.plot(timeNe, Nproj_AU, 'k-', label=r'$N_e^{proj}(t)$')
  plt.plot(unitx*timePost, Npost_AU, 'b-', label=r'$\int_{\Omega} \left( n_{td}(r,t) - n_{gs}(r) \right)^{+} d^{3}r$')
#  plt.legend()

for xc in xposition:
  plt.axvline(x=xc, color='k', linestyle='--')
# plt.plot(unitx*timePost, 10.*Npost_2in_AU, 'r-', label=r'$\int_{\Omega} \left( n_{td}(r,t)) - n_{gs}(r) \right)^2 d^{r}$')
# plt.plot(unitx*timePost, 1./80. * Npost_2out_AU, 'g-', label=r'$\left[ \int_{\Omega} ( N_{td}(r,t) - N_{gs}(r) ) d^{3}r \right]^{2}$')
#print timeCurrent.shape, Nnew.shape
# plt.plot(timeCurrent, (Nnew/10.26**3 * dtCurrent), 'g-', label=r'$\int_{-\infty}^{t} J(t^{\prime}) A(t^\prime) dt^{\prime} $') #TODO: still a problem in units
plt.xlabel('Time (fs)')
# plt.ylabel('Density (a.u.)')
# plt.legend(loc=2)
plt.legend(bbox_to_anchor=(0, -0.3, 1, 0), loc='right', ncol=2, mode="expand", borderaxespad=0)
plt.savefig('KeldyshAnalysis.eps', bbox_inches='tight')
plt.savefig('KeldyshAnalysis.png', bbox_inches='tight')
# plt.show()

print "** Info: Script ended sucessfully."
print "** Goodbye."
