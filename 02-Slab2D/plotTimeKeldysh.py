#!/usr/bin/env python
#-*- coding: utf-8 -*-

# @module plotTimeKeldysh
# Plots the TD-n_excited density of electrons, automatically compares with Keldysh model, assuming a limiter. 
# Plots the TD-energy. 
# Plots the TD-fields and currents. 

## IMPORT LIBRARIES
import numpy as np
from numpy import genfromtxt, loadtxt, chararray, fft
# from pylab import *
import sys
import decimal

import calculateVolume as Sample
import libKeldysh as Keldysh
from libAtomicUnits import *
from libPlottingSettings import *
from libString import *


## TODO: separate plotting and density post-calculation into libs

## Plotting parameters
plotA=311; plotB=312; plotC=313
# Size of main canvas
SizeX = 10
SizeY = 15
RatioA = 10 #(decreases A amplitude visually by RatioA times.

## User choices
UnitSystem="SI" #Defines the system unit for the observables. Time will remain in SI. 
QuantityOfInterest = "density" #number" #Would you prefer to work densities (density), or particle number (number)?
EnableLogScale = 0
ShowFinalPlot  = 1

## Load geometry, n_el, laser and total_current files (atomic units)
geometry = "../Si.xyz" #static method
Ne=loadtxt("../td.general/n_ex", dtype='str', skiprows=9) #total number of excited electrons, not a density!
#laser=loadtxt("../td.general/laser", dtype='str', skiprows=6) #should be replaced by genfromtxt method, as some errors are present in laser file
laser=np.genfromtxt("../td.general/laser", dtype='str', skip_header=9, loose=False, invalid_raise=True, filling_values=0.)
# PostTreatedDensity = loadtxt("postCalculationNexcited.log", dtype='str', skiprows=1) #Unit: number of electrons, which should be multiplied by cell volume (AU). 
TotalCurrent=loadtxt("../td.general/total_current", dtype='str', skiprows=4)
Energy = np.genfromtxt("../td.general/energy", dtype='str', skip_header=5)

## Get the vector potential from laser file
### Building laser profile...
time = np.asfarray(laser[:,1]) #column 2
A1x = np.asfarray(laser[:,2]) #column number 3
A2x = np.asfarray(laser[:,5]) #column number 6 

# c_AU = 137.035999139 #TODO: transfer to lib and remove this line

timestepTable = np.diff(time) #size N-1
TimeStep = timestepTable[0] #they are all equal, just select the first
A_squared = A1x**2
A_integrated = TimeStep * np.cumsum(A1x) #int_{-\inf}^{t} A(t') dt'

Efield = - (np.diff(A1x)/np.diff(time)) #correct value is *1/c 
Efield2= - (np.diff(A2x)/np.diff(time))

Efield_amplitude_SI=np.max(abs(Field_AU_to_SI(Efield / c_AU)))
print "** Info: Electric field amplitude = "+str(Efield_amplitude_SI / 1E9)+" V/nm."
 
### N_Keldysh and N_Gruzdev
print "Getting input from octopus..."
wavelength_inp = LookUpStringInFile("../inp", "lambdaA =")
wavelength_SI = float(wavelength_inp[10:17])
PulseDuration_inp = LookUpStringInFile("../inp", "tauAsi = ")
PulseDuration_SI = float(PulseDuration_inp[8:15])
DoublePulseDelay_inp = LookUpStringInFile("../inp", "pulseDelaySI = ")
DoublePulseDelay_SI = float(DoublePulseDelay_inp[15:22])

# Explicit data
# PulseDuration_SI = 10E-15
Egap_SI = 2.56*e #[Waroquiers and Gonze, Phys. Rev. B 87, 075121 (2013)]

meff = 0.2226 #[Laflamme Janssen and X. Gonze, Physical Review B 93, 205147 (2016)]
#0.18 #Sokolowski-Tinten et al, Phys. Rev. B 61, 2643 (2000)

##==== KELDYSH CALCULATION ====
ValenceBandMax = 5E28 
dtPrecision = 1E-17
order=50
CentralTime = 2.*PulseDuration_SI #20E-15
ShowPlot=False
PeakIntensity_SI = Keldysh.FieldToIntensity(Efield_amplitude_SI)
PeakFluence_SI = PeakIntensity_SI * PulseDuration_SI
PeakFluenceText=round(PeakFluence_SI*1E-4*1E3)
PeakIntensityText = PeakIntensity_SI * 1E-4 #W/cm2
print "** Info: equiv. Peak fluence = "+str(PeakFluenceText)+" mJ/cm2."
print "** Info: equiv. Peak intensity = "+str(PeakIntensityText)+" W/cm2."
timeKeldysh, N_Keldysh_SI, N_Gruzdev_SI = Keldysh.plotPulseToDensity(Egap_SI, meff, wavelength_SI, PulseDuration_SI, PeakFluence_SI, dtPrecision, order, ShowPlot, CentralTime, ValenceBandMax)  

##==== Prepare plotting of N_EX (projection method)  ====
timeNe = np.asfarray(Ne[:,1])
Nproj_AU = np.asfarray(Ne[:,2]) #number of particles, not a density!

##==== Prepare plotting of energies ====
timeEnergies   = np.asfarray(Energy[:,1]) #atomic units (! 1 is 2nd column)
ElectronEnergy = np.asfarray(Energy[:,5]) #Hartree
TotalEnergy    = np.asfarray(Energy[:,2]) #Hartree
IonEnergy      = np.asfarray(Energy[:,3]) #Hartree
Hartree        = np.asfarray(Energy[:,7]) #Hartree

if(timeEnergies.size != ElectronEnergy.size):
  print "** Error: sizes for energies are different. "
  print timeEnergies.size, ElectronEnergy.size
  exit()

# Normalizing to 0
ElectronEnergy = ElectronEnergy - ElectronEnergy.min()
TotalEnergy    = TotalEnergy - TotalEnergy.min()
IonEnergy      = IonEnergy - IonEnergy.min()
Hartree        = Hartree - Hartree.min()

##===== Prepare plotting of total_current ======
ShrinkBy=3

timeCurrent=np.asfarray(TotalCurrent[:,1])
# timeCurrent=timeCurrent[:timeCurrent.size-ShrinkBy]
maxSize=timeCurrent.size-ShrinkBy

CurrentX=np.asfarray(TotalCurrent[:,2])
# CurrentX=CurrentX[:maxSize]

print "Current: "+str(timeCurrent.size)+", "+str(CurrentX.size)

Efield_AU =Efield/c_AU #fixing dimension #TODO: clean me that...
Efield2_AU=Efield2/c_AU

##== Fourier transforms ==
# Efield_AU_FT  =fft.rfft(Efield[:maxSize]/c_AU)
# CurrentX_AU_FT=fft.rfft(CurrentX[:maxSize])
# Conductivity_Freqs=fft.rfft(timeCurrent[:maxSize])

##TODO: use SI units before transformations ?
# print "** Conductivity"
# Epsilon_maxSize = CurrentX_AU_FT.size
# # print Efield_AU_FT.size
# Conductivity = np.divide(CurrentX_AU_FT, Efield_AU_FT) #working in AU space ...
# Conductivity_plotX = Conductivity_Freqs #Conductivity_Freqs # timeCurrent[1:maxSize]
# Conductivity_plotY = Efield_AU_FT # Conductivity # Efield_AU[1:maxSize]
# print "Conductivity_X: "+str(Conductivity_plotX.size)
# print "Conductivity_Y: "+str(Conductivity_plotY.size)


##=== POST-PROCESSED charge density (n_ex from n_tot) ===
# timePost = np.asfarray(PostTreatedDensity[:,2]) #3: time SI
# NePost = np.asfarray(PostTreatedDensity[:,3]) #4: N_e number of electrons to multiply by cell volume (AU)
# NePostSquaredInside = np.asfarray(PostTreatedDensity[:,6]) #7: N_e squared inside
##  NePostSquaredOutside = NePost**2 #N_e squared outside the integral

##=== CALCULATION OF THE Single + Supercell Volume (at. u.)
SuperCellVolume_AU = Sample.VolumeOfSupercell_Log("../")
SuperCellVolume_SI = Volume_AU_to_SI(SuperCellVolume_AU)
SampleCellVolume_AU = SuperCellVolume_AU - Sample.VolumeOfVacuum("../") #Volume of the sample (without vacuum region) 
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

##====== Excited density spatial integration ========

# Calculation using quantum projections < phi_gs | phi(t) >
#Nproj_AU = density #* CellVolume_AU 
# Nproj_AU = Nproj_AU * CellVolume_AU #Should we really multiply by volume ?
density_proj_AU = Nproj_AU / SampleCellVolume_AU #This is now a density (a.u.)
# Calculation of particle number from charge density (octopus.output.density)
# Npost_AU = NePost * CellVolume_AU #NOW, it's a number of electrons. :)
# Npost_2in_AU = NePostSquaredInside * CellVolume_AU #Same
# Npost_2out_AU= NePostSquaredOutside * CellVolume_AU #Same

# density_post_AU = Npost_AU / SuperCellVolume_AU

# print CurrentX.size, A1x.size
# Shrink A1x to size of CurrentX


##===== USED ? Developed when attempting to compute N_ex(N_tot) ==== 
### Attempt using J(r,t) x A(r,t). 
Overlap=np.array([CurrentX.size, A1x.size])
indexmax=Overlap.min()
A1xShrinked = A1x[0:indexmax]
A2xShrinked = A2x[0:indexmax]
CurrentXShrinked=A1x[0:indexmax]
NnewDerivDt=np.multiply(CurrentXShrinked,A1xShrinked) #still to multiply by dt and sum.
Nnew = np.cumsum(NnewDerivDt)

##====== Converting timescales from AU to SI ============
UnitTime = 1E15 #we show time in fs
TimeAUtoSI = Time_AU_to_SI(1.)
time = UnitTime * time * TimeAUtoSI
timeNe = UnitTime * timeNe * TimeAUtoSI

timeCurrentdt = timeCurrent #* TimeAUtoSI
timeCurrent=UnitTime * timeCurrent * TimeAUtoSI
dtCurrent = timeCurrentdt[1]-timeCurrentdt[0]

##=== Convert density from AU to SI =====
# Avoid ill-defined units
if(UnitSystem == "SI"):
  print "** Info: SI system chosed."
  observable_proj = Density_AU_to_SI(density_proj_AU)
#  observable_post = Density_AU_to_SI(density_post_AU)
elif(UnitSystem == "AU"):
  print "** Info: AU system chosed."
  observable_proj = density_proj_AU
#  observable_post = density_post_AU

else: 
  print "** Error: Unit system was ill defined. Set UnitSystem=[SI,AU] in plotTimeKeldysh.py."
  exit()



### Converting fields from AU to SI

print "TimeStep = "+str(dtCurrent)

## ================================== TIME FOR PLOTTING ====================
### Defining plotting boundaries
xmin=0.; xmax=50.; unitx=1E15
ymax=2E28; #max for density plot
ymax2=1E0; #max for field plot

## Defining plot structure
plt.figure(figsize=(SizeX,SizeY))
ax1 = plt.subplot(plotA)

## Plotting electric field
plt.title(r"$\lambda=$"+str(wavelength_SI*1E9)+" nm, $I_0=$"+str('%.2E' % decimal.Decimal(PeakIntensityText))+" W/cm$^2$, delay="+str(float(DoublePulseDelay_SI*1E15))+" fs")
#plt.xlabel('Time (fs)')
# plt.ylabel('Vector potential (a.u.)')
plt.ylabel('Field amplitude (a.u.)')
# plt.xlim((0, time.max()))
plt.xlim((xmin,xmax))
# plt.ylim((-ymax2,ymax2))

# plt.plot(time[1:], c*EtotX, 'k-', label='c x E(t)')
plt.plot(time, A1x/RatioA, 'b-', label=r'$A_1(t)$/'+str(RatioA)) #GOOD
plt.plot(time, A2x/RatioA, 'b-', label=r'$A_2(t)$/'+str(RatioA))
# plt.plot(time, A_integrated, 'r-', label="$\int_{-\infty}^{t} A(t') dt'$")
# plt.plot(time, A_squared, 'g-', label="A^2 (t)")
plt.plot(time[1:], Efield_AU*c_AU, 'k-', label=r"$c E_1(t)$")
plt.plot(time[1:], Efield2_AU*c_AU, 'k-', label=r"$c E_2(t)$")

# plt.plot(timeCurrent, TimeStep * np.cumsum(CurrentX), 'r--', label=r'$P(t)$')
plt.plot(timeCurrent, CurrentX, 'r-', label=r'$J(t)$')

# wavelength_SI is now automatic
wavelength_AU = length_SI_to_AU(wavelength_SI)

##Let's add a grid to compare the A extrema and Ne(t) extrema
m=np.arange(0,200) #[-5.,-4.,-3.,-2.,-1.,0.,1.,2.,3.,4.,5.]
xposition = np.multiply(m,unitx*wavelength_SI/(2.*c))-3.*pi/2.
for xc in xposition:
  plt.axvline(x=xc, color='k', linestyle=':')
plt.legend(loc=4)

### Plotting density evolution
ax2 = plt.subplot(plotB)
# plt.xlim((0,unitx*timePost.max()))
plt.xlim((xmin,xmax))
# plt.ylim((0,ymax))

if(QuantityOfInterest == "density"):
  if(UnitSystem=="AU"):
    plt.ylabel('Density exc. el. (a.u.)')
  elif(UnitSystem=="SI"):
    plt.ylabel(r'Density exc. el. (m$^{-3}$)')
  plt.plot(timeNe, observable_proj, 'k-', label=r'$n_e^{proj}(t)$')
#  plt.plot(unitx*timePost, observable_post, 'k--', label=r'TD-DFT: \frac{1}{V_{tot}} $\int_{\Omega} \left( n_{td}(r,t) - n_{gs}(r) \right)^{+} d^{3}r$')
  plt.plot(timeKeldysh*unitx, N_Keldysh_SI, 'r-', label=r'[Keldysh (1964)]')
  plt.plot(timeKeldysh*unitx, N_Gruzdev_SI, 'b-', label=r'[Gruzdev (2014)]')
elif (QuantityOfInterest == "number"):
  plt.ylabel('Number of exc. electrons')
  plt.plot(timeNe, Nproj_AU, 'k-', label=r'$N_e^{proj}(t)$')
#  plt.plot(unitx*timePost, Npost_AU, 'b-', label=r'$\int_{\Omega} \left( n_{td}(r,t) - n_{gs}(r) \right)^{+} d^{3}r$')
#  plt.legend()

## Plotting vertical bars for comparison A(t) and Ne(t)
for xc in xposition:
  plt.axvline(x=xc, color='k', linestyle=':')
# plt.plot(unitx*timePost, 10.*Npost_2in_AU, 'r-', label=r'$\int_{\Omega} \left( n_{td}(r,t)) - n_{gs}(r) \right)^2 d^{r}$')
# plt.plot(unitx*timePost, 1./80. * Npost_2out_AU, 'g-', label=r'$\left[ \int_{\Omega} ( N_{td}(r,t) - N_{gs}(r) ) d^{3}r \right]^{2}$')

# plt.plot(timeCurrent, (Nnew/10.26**3 * dtCurrent), 'g-', label=r'$\int_{-\infty}^{t} J(t^{\prime}) A(t^\prime) dt^{\prime} $') #TODO: still a problem in units
plt.xlabel('Time (fs)')
# plt.ylabel('Density (a.u.)')
# plt.legend(loc=2)
# plt.legend(bbox_to_anchor=(0, -0.3, 1, 0), loc='right', ncol=2, mode="expand", borderaxespad=0)

## ============ PLOTTING ENERGIES ======
ax3 = plt.subplot(plotC)
plt.xlim((xmin,xmax))
plt.plot(UnitTime * TimeAUtoSI * timeEnergies, ElectronEnergy, 'k-', label=r'$\xi_{el} (tot)$')
plt.plot(UnitTime * TimeAUtoSI * timeEnergies, IonEnergy,      'b--', label=r"$\xi_{ion}$")
# plt.plot(UnitTime * TimeAUtoSI * timeEnergies, Hartree,      'r-', label=r"$\xi_{Hartree}$")
plt.xlabel('Time (fs)')
plt.ylabel('Energy [Ha]')
plt.legend(loc=2)
## Plotting vertical bars for comparison A(t) and Ne(t)
for xc in xposition:
  plt.axvline(x=xc, color='k', linestyle=':')

## ============ PLOTTING FOURIER SPECTRA =========
# ax4 = plt.subplot(414)
# plt.plot(Conductivity_plotX, Conductivity_plotY, 'ko', label=r'$\sigma(\omega)$')
# ax3.set_xscale("log", nonposx='clip')

## ======== Applying logscales...

if(EnableLogScale == 1):
  ax2.set_yscale("log", nonposy='clip')
  ax2.legend(loc=4)
  ax3.set_yscale("log", nonposy='clip')
  ax3.legend(loc=4)
else:
  ax2.legend(loc=4)

## ==== Writing final picture =====

filename='KeldyshAnalysis'
if(EnableLogScale==1): 
  filename=filename+"-LogScale"
plt.savefig(filename+'.eps', bbox_inches='tight')
plt.savefig(filename+'.png', bbox_inches='tight')
if(ShowFinalPlot==1):
  plt.show()

print "** Info: Script ended sucessfully."
print "** Goodbye."
