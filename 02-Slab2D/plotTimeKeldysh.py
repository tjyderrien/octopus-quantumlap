#!/usr/bin/env python2
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2018 T. J.-Y. Derrien
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
from libPlotting import *


## Materials properties
Egap_SI = 2.56*e #[Waroquiers and Gonze, Phys. Rev. B 87, 075121 (2013)] #TODO: this will depend on which functional we have used! TB09 or LDA ! 
meff = 0.2226 #[Laflamme Janssen and X. Gonze, Physical Review B 93, 205147 (2016)] #SAME !
#0.18 #Sokolowski-Tinten et al, Phys. Rev. B 61, 2643 (2000)

## Plotting: parameters for controlling the canvas
## 3 canvas
#plotA=311; plotB=312; plotC=313;
## 4 canvas
plotA=411; plotB=412; plotC=413; plotD=414
# Size of main canvas
SizeX = 10
SizeY = 20
RatioA = 10 #(decreases A amplitude visually by RatioA times.

## Position of the captions
CaptionX = 1.00
CaptionY = 1.00

## User choices
UnitSystem="SI" #Defines the system unit for the observables. Time will remain in SI. 
QuantityOfInterest = "density" #number" #Would you prefer to work densities (density), or particle number (number)?
EnableLogScale = 2 #1: yes, 2: auto from command line

ShowFinalPlotDefault  = 0 #Forcing the value of ShowFinalPlot
PlottingMode          = "auto" # "doublepulse" #nonthermalmelting, keldysh #TODO: define automatically as function of parameters
PlottingMethod = "detect" #bruteforce

maxHHG = 3*200 #maximum number of samples to cut the spectrum with

Header="[octopus-slabs] "

if(PlottingMode == "auto"):
## METHOD BRUTEFORCE
  if(PlottingMethod == "bruteforce"): #{{{
    PlottingMode = sys.argv[1]
    if( not (PlottingMode == "nonthermalmelting" or PlottingMode == "doublepulse" or PlottingMode == "keldysh") ):
      print "Usage: plotTimeKeldysh.sh <PlottingMode = nonthermalmelting | doublepulse | keldysh> <LogScale = 0 | 1> <ShowFinalPlot = 0 | 1>"
      exit()
  #}}}
## METHOD detection in input file
  else: #{{{
    filename = "../inp"
    Efield2_inp = LookUpStringInFile(filename, "EfieldBsi =")
    Efield2_SI  = float(Efield2_inp[11:15])
    MoveIonChain = LookUpStringInFile(filename, "MoveIons = ")
    MoveIonValue = MoveIonChain[11:15]
    MoveIonClean = MoveIonValue.strip()
    LastNexcLine = GetLastLine("../td.general/n_ex") 
    LastNexcLineArray = np.array([ float(x) for x in LastNexcLine.split() ])
    
    if(Efield2_SI > 1E0): #double pulse enabled
      PlottingMode="doublepulse" #priority 1
    # elif(Time_AU_to_SI(LastNexcLineArray[1]) > 100e-15 ): #latest time of simulation is longer than 100 fs
#TODO: nonthermalmelting is meaningful only if MoveIon=yes... 
    elif(MoveIonClean == "yes"):
      PlottingMode="nonthermalmelting" #priority 2
    else: 
      PlottingMode="keldysh" #priority 3
    print "** Plotting mode automatic: "+PlottingMode
  #}}}
if(int(EnableLogScale) == 2):
  try: 
    EnableLogScale = int(sys.argv[2]) #read from terminal
#    ShowFinalPlot  = 0     # MORE SAFE. 
  except: 
    print "Usage: plotTimeKeldysh.py <PlottingMode = auto | keldysh | doublepulse | nonthermalmelting> <LogScale? 0|1> <ShowFinalPlot? 0 | 1>."
    exit()
if(sys.argv[3] == ""):
  ShowFinalPlot = 0
elif(sys.argv[3] == "1"):
  ShowFinalPlot = 1
else:
  ShowFinalPlot = ShowFinalPlotDefault
  

geometry, Ne, laser, TotalCurrent, Energy = LoadRunOutputs("")

# Treat the data generated with Octopus. 

## Converts a vector potential A to electric field 
#def GaugeFieldToElectricField(laser):

## Get the vector potential from laser file
### Building laser profile...
time = np.asfarray(laser[:,1]) #column 2
A1x  = np.asfarray(laser[:,2]) #column number 3
A1y  = np.asfarray(laser[:,3]) #column number 4
A2x  = np.asfarray(laser[:,5]) #column number 6 
A2y  = np.asfarray(laser[:,6]) #column number 7

timestepTable = np.diff(time) #size N-1
TimeStep = timestepTable[0] #they are all equal, just select the first
#A_squared = A1x**2
#A_integrated = TimeStep * np.cumsum(A1x) #int_{-\inf}^{t} A(t') dt'

# oscillating fields from octopus (real valued)
Efield_X  = - (np.diff(A1x)/np.diff(time)) #physical value is *1/c 
Efield2_X = - (np.diff(A2x)/np.diff(time)) #physical value is *1/c
Efield_Y  = - (np.diff(A1y)/np.diff(time)) #physical value is *1/c 
Efield2_Y = - (np.diff(A2y)/np.diff(time)) #physical value is *1/c

## Correcting the value by 1/c_AU
Efield_X_AU =Efield_X/c_AU
Efield_Y_AU =Efield_Y/c_AU
Efield2_X_AU=Efield2_X/c_AU
Efield2_Y_AU=Efield2_Y/c_AU
Efield_X_SI  = Field_AU_to_SI(Efield_X_AU)
Efield2_X_SI = Field_AU_to_SI(Efield2_X_AU)

# E field amplitudes (in atomic units, positive real-valued)
Efield_norm_AU  = np.sqrt(Efield_X_AU**2  + Efield_Y_AU**2) 
Efield2_norm_AU = np.sqrt(Efield2_X_AU**2 + Efield2_Y_AU**2)
Efield_norm_SI  = Field_AU_to_SI(Efield_norm_AU)
Efield2_norm_SI = Field_AU_to_SI(Efield2_norm_AU)
Efield_amplitude_SI =np.max(abs(Efield_norm_SI))
Efield_amplitudeB_SI=np.max(abs(Efield2_norm_SI))
print Header+"** Info: Electric field 1: amplitude = "+str(Efield_amplitude_SI / 1E9)+" V/nm."
print Header+"** Info: Electric field 2: amplitude = "+str(Efield_amplitudeB_SI / 1E9)+" V/nm."

### Computing the total oscillating field for mixed pulses (from octopus, real-valued)
EfieldAB_X_AU = Efield_X_AU + Efield2_X_AU
EfieldAB_Y_AU = Efield_Y_AU + Efield2_Y_AU
EfieldAB_amplitude_AU = np.sqrt(EfieldAB_X_AU**2 + EfieldAB_Y_AU**2) #if polarization are crossed, Keldysh will NOT WORK ! But this is the point: generating the data we need in that case :) 
EfieldAB_amplitude_SI = np.max(abs(Field_AU_to_SI(EfieldAB_amplitude_AU)))
print Header+"** Info: Electric field 1+2: amplitude = "+str(EfieldAB_amplitude_SI / 1E9)+" V/nm."
 
### Importing input from octopus 
print Header+"== Getting input from octopus... =="
print Header+" PULSE 1 ..."
wavelength_inp       = LookUpStringInFile("../inp", "lambdaA =")
wavelength_SI        = float(wavelength_inp[10:17])
PulseDuration_inp    = LookUpStringInFile("../inp", "tauAsi = ")
PulseDuration_SI     = float(PulseDuration_inp[8:15])

print Header+" PULSE 2 ..."
wavelengthB_inp      = LookUpStringInFile("../inp", "lambdaB =")
wavelengthB_SI       = float(wavelengthB_inp[10:17])
PulseDurationB_inp   = LookUpStringInFile("../inp", "tauBsi = ")
PulseDurationB_SI    = float(PulseDurationB_inp[8:15])
DoublePulseDelay_inp = LookUpStringInFile("../inp", "pulseDelaySI = ")
DoublePulseDelay_SI  = float(DoublePulseDelay_inp[15:22])

# Gaussian | SinSquared pulse | Soft Top Hat pulse?
PulseType_inp   = LookUpStringInFile("../inp", "FieldColorA-")
PulseType_value = PulseType_inp[49:] #this allows to have unique solution
if("FieldColorA-Gaussian" in PulseType_value):
  print Header+"** Info: Gaussian pulse detected..."
  PulseType="gaussian"
elif ("FieldColorA-SinSquared" in PulseType_value):
  print Header+"** Info: Squared Sinus pulse detected..."
  PulseType="SinSquared"
elif ("FieldColorA-STH" in PulseType_value):
  print Header+"** Info: Soft Top Hat pulse detected..."
  PulseType="STH"
elif ("FieldColorA-StH" in PulseType_value):
  print Header+"** Info: Soft Top Hat pulse [corrected] detected..."
  PulseType="StH"
else: 
  print Header+"** Error: pulse type was not recognized."
  print Header+"** Error value: PulseType_value = "+str(PulseType_value)
  exit()

##==== KELDYSH CALCULATION PULSE 1 ==== NOTE: /!\ ONLY VALID FOR GAUSSIAN PULSES !!
ValenceBandMax = 5E28 
dtPrecision = 1E-17
order=50
CentralTime = 2.*PulseDuration_SI #20E-15
ShowPlot=False
PeakIntensity_SI  = Keldysh.FieldToIntensity(Efield_amplitude_SI)
PeakFluence_SI    = PeakIntensity_SI * PulseDuration_SI
PeakFluenceText   = round(PeakFluence_SI*1E-4*1E3)
PeakIntensityText = PeakIntensity_SI * 1E-4 #W/cm2
print Header+"** Info: equiv. Peak fluence = "+str(PeakFluenceText)+" mJ/cm2."
print Header+"** Info: equiv. Peak intensity = "+str(PeakIntensityText)+" W/cm2."

##==== KELDYSH CALCULATION PULSE 2 ==== NOTE: /!\ ONLY VALID FOR GAUSSIAN PULSES !!
PeakIntensityB_SI = Keldysh.FieldToIntensity(np.max(Efield2_norm_SI))
PeakFluenceB_SI = PeakIntensityB_SI * PulseDurationB_SI
PeakFluenceBText=round(PeakFluenceB_SI*1E-4*1E3)
PeakIntensityBText = PeakIntensityB_SI * 1E-4 #W/cm2
print Header+"** Info: equiv. Peak fluence = "+str(PeakFluenceBText)+" mJ/cm2."
print Header+"** Info: equiv. Peak intensity = "+str(PeakIntensityBText)+" W/cm2."

CentralTimeB = 2.*PulseDurationB_SI + DoublePulseDelay_SI #20E-15
# if(PlottingMode == "doublepulse"):

##==== Building the laser pulses for Keldysh ======
##(they cannot be the ones of octopus, as we need
## the complex space to build the bicolor envelope)
if(PulseType == "SinSquared"): 
  t0   =  PulseDuration_SI #initial time is actually 0 in octopus. 
  tmin = -PulseDuration_SI+t0 #zero...
  if(PlottingMode=="doublepulse"):
    tmax =  PulseDuration_SI + DoublePulseDelay_SI + PulseDurationB_SI + t0
  else:
    tmax = PulseDuration_SI + t0
elif(PulseType == "gaussian"): 
  t0   =  2.*PulseDuration_SI
  tmin = -2.*PulseDuration_SI + t0
  tmax =  2.*PulseDuration_SI + t0
elif(PulseType == "STH" or PulseType=="StH"): 
  t0   =  2.*PulseDuration_SI
  tmin = -2.*PulseDuration_SI + t0
  tmax =  2.*PulseDuration_SI + t0 #BUG
else: 
  print Header+"** Error: pulse type was not recognized. "
  print Header+"** Error details: PulseType = "+str(PulseType)
  exit()
  
instants=np.arange(tmin,tmax,dtPrecision)

CEP1=0.; CEP2=0. #TODO: introduce CEP1 and CEP2 auto reading from octopus input file

print Header+"Defining Pulses for Keldysh analysis..."
if(PulseType == "SinSquared"): 
  FieldEnvelope,  FieldComplex    = Keldysh.PulseSquaredSinTemporalShape(instants, PulseDuration_SI,  Efield_amplitude_SI,  wavelength_SI, CEP1, t0, 0.) 
  FieldEnvelope2, FieldComplex2   = Keldysh.PulseSquaredSinTemporalShape(instants, PulseDurationB_SI, Efield_amplitudeB_SI, wavelength_SI, CEP2, t0, DoublePulseDelay_SI) 
  FieldEnvelopeAB, FieldComplexAB = Keldysh.PulseSquaredSinTemporalShapeDoublePulse(instants, PulseDuration_SI, PulseDurationB_SI, Efield_amplitude_SI, Efield_amplitudeB_SI, wavelength_SI, wavelengthB_SI, CEP1, CEP2, t0, DoublePulseDelay_SI) 
elif(PulseType == "gaussian"):
  IntensityEnvelope  = Keldysh.PulseGaussianTemporalShape(instants, PulseDuration_SI,  Keldysh.FieldToIntensity(Efield_amplitude_SI), t0)
  FieldEnvelope      = Keldysh.IntensityToField(IntensityEnvelope) 
  IntensityEnvelope2 = Keldysh.PulseGaussianTemporalShape(instants, PulseDurationB_SI, Keldysh.FieldToIntensity(Efield_amplitudeB_SI), t0)
  FieldEnvelope2     = Keldysh.IntensityToField(IntensityEnvelope2)
  #TODO: add double-gaussian-pulses here
elif(PulseType == "STH" or PulseType=="StH"): 
  #here, we don't compute the density but just output the W_PI, because there is not analytical solution for dNe/dt. 
  #As a result, we only put one instant into the pulse to obtain a scalar at the peak intensity. 
  #IntensityEnvelope  = Keldysh.PulseGaussianTemporalShape(t0, PulseDuration_SI,  Keldysh.FieldToIntensity(Efield_amplitude_SI), t0) #BUG: we just need the peak amplitude of the field
  FieldEnvelope      = Efield_amplitude_SI #Keldysh.IntensityToField(IntensityEnvelope)
  #IntensityEnvelope2 = Keldysh.PulseGaussianTemporalShape(t0, PulseDurationB_SI, Keldysh.FieldToIntensity(Efield_amplitudeB_SI), t0)
  FieldEnvelope2     = Efield2_norm_SI #Keldysh.IntensityToField(IntensityEnvelope2)
else: 
  print Header+"** Error in definition of the pulses for Keldysh. Pulse was not recognized."
  exit()
  
if(PulseType == "gaussian" or PulseType=="SinSquared"): #{{{
# In that case, we do have analytical solutions for 1 color. 
  ## Computing final excited electron density for Pulse 1
  print Header+"Computing Keldysh for Pulse 1: "
  timeKeldysh, N_Keldysh_SI, N_Gruzdev_SI = Keldysh.plotPulseToDensity(Egap_SI, meff, wavelength_SI, PulseDuration_SI, FieldEnvelope, dtPrecision, order, ShowPlot, CentralTime, ValenceBandMax)  
  PeakWpi_Gruzdev = Keldysh.GenerateKeldyshDatabase(Egap_SI, meff, wavelength_SI, FieldEnvelope, order) #NOTE: is this really useful? 
  ## Computing final excited electron density for Pulse 2
  if(PlottingMode=="doublepulse"):
    print Header+"Computing Keldysh for Pulse 2: "
    try: 
      timeKeldyshB, N_KeldyshB_SI, N_GruzdevB_SI = Keldysh.plotPulseToDensity(Egap_SI, meff, wavelengthB_SI, PulseDurationB_SI, FieldEnvelope2, dtPrecision, order, ShowPlot, CentralTimeB, ValenceBandMax)
    except: 
      N_KeldyshB_SI = N_Keldysh_SI * 0e0
      N_GruzdevB_SI = N_Gruzdev_SI * 0e0

    ##==== KELDYSH CALCULATION FOR BOTH PULSES =====
    print Header+"Preparing double pulse for Keldysh..."
  #   PeakIntensityAB_SI  = Keldysh.FieldToIntensity(EfieldAB_amplitude_SI)
  #   PeakIntensityABText = PeakIntensityAB_SI * 1E-4 #W/cm2

    if(wavelength_SI == wavelengthB_SI):
  #    IntensityEnvelop = PeakIntensityAB_SI # Keldysh.FieldToIntensity(EfieldAB_amplitude_SI)
      timeKeldyshAB, N_KeldyshAB_SI, N_GruzdevAB_SI = Keldysh.plotPulseToDensity(Egap_SI, meff, wavelength_SI, PulseDurationB_SI+CentralTimeB, FieldEnvelopeAB.real, dtPrecision, order, ShowPlot, CentralTimeB, ValenceBandMax)
    else:
      print "Bicolor Keldysh is not implemented! See V. Zhukov [Novosibirsk University]"
#}}}
elif (PulseType == "STH" or PulseType == "StH"): 
  # For this case, we only output the w_pi. 
  PeakWpi_Gruzdev = Keldysh.GenerateKeldyshDatabase(Egap_SI, meff, wavelength_SI, FieldEnvelope, order)
else: 
  print Header+"** Error: computation of density | excitation rate did not recognize the shape of the pulse."
  exit()

##==== Prepare plotting of N_EX (projection method)  ====
timeNe = np.asfarray(Ne[:,1])
Nproj_AU = np.asfarray(Ne[:,2]) #number of particles, not a density!

##==== Prepare plotting of energies ====
timeEnergies      = np.asfarray(Energy[:,1]) #atomic units (! 1 is 2nd column)
TotalEnergy       = np.asfarray(Energy[:,2]) #Hartree
IonEnergy         = np.asfarray(Energy[:,3]) #Hartree, Kinetic energy !
#NOTE: We dont need the ion-ion when matter does not move. 
#TODO: But we DO need it when MoveIon = yes! 
ElectronEnergy    = np.asfarray(Energy[:,5]) #Hartree
EigenValuesEnergy = np.asfarray(Energy[:,6]) #Hartree
Hartree           = np.asfarray(Energy[:,7]) #Hartree energy (electron-electron interaction)
IntVxcNEnergy     = np.asfarray(Energy[:,8]) #Hartree, int ( vxc n(r) )
ExchangeEnergy    = np.asfarray(Energy[:,9]) #Hartree, exchange 
CorrelationEnergy = np.asfarray(Energy[:,10]) #Hartree, correlation
try:
  ElectronKineticEnergy = np.asfarray(Energy[:,11]) #Hartree, electron kinetic energy
except:
  ElectronKineticEnergy = 0. #NOTE: was not outputed before 2017 Nov 29th, and might be even wrong. 

# Somehow checking the data integrity
if(timeEnergies.size != ElectronEnergy.size):
  print "** Error: sizes for energies are different. "
  print timeEnergies.size, ElectronEnergy.size
  exit()

# Controlling that we well understood the energy balance. 
# RebuildTotalEnergy = EigenValuesEnergy - Hartree + ExchangeEnergy + CorrelationEnergy - IntVxcEnergy
# if(RebuildTotalEnergy != TotalEnergy):
#   print Header+"** Warning: Kinetic energy may not be well reconstructed. "

# Normalizing to 0
ElectronEnergy = ElectronEnergy - ElectronEnergy.min()
TotalEnergy    = TotalEnergy - TotalEnergy.min()
IonEnergy      = IonEnergy - IonEnergy.min()
Hartree        = Hartree - Hartree.min()

##===== Prepare plotting of total_current ======
ShrinkBy=0

timeCurrent=np.asfarray(TotalCurrent[:,1])
# timeCurrent=timeCurrent[:timeCurrent.size-ShrinkBy]
maxSize=timeCurrent.size-ShrinkBy

CurrentX=np.asfarray(TotalCurrent[:,2])
# CurrentX=CurrentX[:maxSize]

print "Current: "+str(timeCurrent.size)+", "+str(CurrentX.size)

##== Fourier transforms ==
# Efield_AU_FT  =fft.fft(Efield[:maxSize]/c_AU)
# CurrentX_AU_FT=fft.fft(CurrentX[:maxSize])
# DiffCurrentX_AU=np.diff(CurrentX) #TODO: NTD: do that in Fourier space (much more precise!)
HHG_AU=fft.fft(CurrentX)
HHG_AU_omega=2.*pi*fft.fftfreq(maxSize,d=TimeStep*Time_AU_to_SI(1.))
# c_SI = Velocity_AU_to_SI(c_AU)
HHG_AU_order=HHG_AU_omega/(2.*pi*c/wavelength_SI)
HHG_AU_duration=1./HHG_AU_omega
# Conductivity_Freqs=fft.fft(timeCurrent[:maxSize])


##TODO: Prepare the conductivity plots, also! 
## Do not use FFT, see reasons in Numerical Recipes. 
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
density_proj_AU = Nproj_AU / SampleCellVolume_AU #This is now a density (at.u.)
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

### Building up an instantaneous W_PI(t)
omega = 2.*pi*c/wavelength_SI
wPI_TDDFT = np.diff(observable_proj)/np.diff(timeNe/UnitTime) #virtual excitations make it unclear
# wPI_TDDFT_average = wPI_TDDFT.cumsum() / timeNe[1:].size #not good as depends on the final time
NumberOfTimeSteps = int(2.*PulseDuration_SI/(TimeStep*TimeAUtoSI))
wPI_TDDFT_average = wPI_TDDFT[:NumberOfTimeSteps-1] / NumberOfTimeSteps 
wPI_TDDFT_cos     = np.exp( -1.j*omega * timeNe[1:]  )
wPI_TDDFT_env     = wPI_TDDFT * wPI_TDDFT_cos 
wPI_TDDFT_approx  = observable_proj / (2. * PulseDuration_SI)

print Header+"** Final value for <w_PI> = "+str(wPI_TDDFT_average[wPI_TDDFT_average.size-1])
print Header+"** Final value for w_PI_approx = "+str(wPI_TDDFT_approx[wPI_TDDFT_approx.size-1])

### Converting fields from AU to SI

print "TimeStep = "+str(dtCurrent)

## ================================== TIME FOR PLOTTING ====================
# Managing stepwise maximum for Y1
ymin=1E16; ymax=ValenceBandMax; #max for density plot
if(observable_proj.max() < ymax / 50.):
  ymax = ymax / 50.
elif(observable_proj.max() < ymax / 5.):
  ymax = ymax / 5.
elif(observable_proj.max() > ValenceBandMax):
  ymax=4.*ymax;
# Managing stepwise maximum for Y2
ymax2=5E0; 

ymax11value = (Efield_X_SI+Efield2_X_SI) #TODO: this is not optimum when not overlapping.  
ymax11valueMax = ymax11value.max()
ymax12valueMax = CurrentX.max()
#if(ymax2valueMax < ymax2 / 100.):
  #ymax2 = ymax2 / 100.
#elif(ymax2valueMax < ymax2 / 10.):
  #ymax2 = ymax2 / 10.
#elif(ymax2valueMax > ymax2):
  #ymax2=5.*ymax2;
# Managing wpi plotting
ymin_wpi = 1E40
ymax_wpi = 1E50


### Defining plotting boundaries
if(PlottingMode=="nonthermalmelting"):
  xmin=0.; xmax=500.; unitx=1E15
#  ymin=1E16; ymax=5E28; #max for density plot
#  if(observable_proj.max() > 5E28):
#    ymax=4.*ymax; 
#  ymax2=5E0; #max for field plot
  PlotGrid=0
elif(PlottingMode=="doublepulse"):
  xmin=0.; xmax=60.; unitx=1E15
#  ymin=1E16; ymax=4.*5E28; 
#  ymax2=5E0;
  PlotGrid=1 
elif(PlottingMode=="keldysh"):
  xmin=0.; unitx=1E15; xmax=PulseDuration_SI*unitx*2.;
#  ymin=1E16; ymax=5E28;
#  ymax2=5E0;
  PlotGrid=1

## Defining plot structure
plt.figure(figsize=(SizeX,SizeY))
ax1 = plt.subplot(plotA)

## Plotting electric field
plt.title(r"$\lambda=$"+str(wavelength_SI*1E9)+" nm, $I_0=$"+str('%.2E' % decimal.Decimal(PeakIntensityText))+" W/cm$^2$, delay="+str(float(DoublePulseDelay_SI*1E15))+" fs")
#plt.xlabel('Time (fs)')
# plt.ylabel('Vector potential (a.u.)')
plt.ylabel('Field amplitude (V/m)')
# plt.xlim((0, time.max()))
plt.xlim((xmin,xmax))
plt.ylim((-ymax11valueMax,ymax11valueMax))

# plt.plot(time[1:], c*EtotX, 'k-', label='c x E(t)')
### Showing vector potential is only useful to observe shift with materials excitation. 
# plt.plot(time, (A1x+A2x)/RatioA, 'b-', label=r'$A(t)$/'+str(RatioA)) #GOOD 
#TODO: warn if polarization are not similar
# plt.plot(time, A2x/RatioA, 'b-')
# plt.plot(time, A_integrated, 'r-', label="$\int_{-\infty}^{t} A(t') dt'$")
# plt.plot(time, A_squared, 'g-', label="A^2 (t)")
plot_field, = ax1.plot(time[1:], (Efield_X_SI+Efield2_X_SI), 'k-', label=r"$E(t)$") #E ~ i A omega / c
# plt.plot(time[1:], Efield2_AU*c_AU, 'k-')

#let's plot current in a second axis
ax12 = ax1.twinx()
plt.xlim((xmin,xmax))
plt.ylim((-ymax12valueMax, ymax12valueMax))
plt.ylabel(r'$J(t)$')
# plt.plot(timeCurrent, TimeStep * np.cumsum(CurrentX), 'r--', label=r'$P(t)$')
plot_current, = ax12.plot(timeCurrent, CurrentX, 'r-', label=r'$J(t)$ (a.u.)')
ax12.yaxis.label.set_color(plot_current.get_color()) #colorizes the label
ax12.spines["right"].set_edgecolor(plot_current.get_color()) #colorizes the axis
ax12.tick_params(axis='y', colors=plot_current.get_color()) #colorizes the tics and numbers

# Building the shared legend for subfig1
plot1   = [plot_field, plot_current]
labels1 = [l.get_label() for l in plot1]
ax1.legend(plot1, labels1, loc=2)

# wavelength_SI is now automatic
wavelength_AU = length_SI_to_AU(wavelength_SI)

if(PlotGrid==1):
  ##Let's add a grid to compare the A extrema and Ne(t) extrema
  m=np.arange(0,500) #[-5.,-4.,-3.,-2.,-1.,0.,1.,2.,3.,4.,5.]
  xposition = np.multiply(m,unitx*wavelength_SI/(2.*c))-3.*pi/2.
  for xc in xposition:
    plt.axvline(x=xc, color='#b3b3b3', linestyle=':')

# plt.legend(bbox_to_anchor=(CaptionX, 0, loc=1, borderaxespad=0.)
#plt.legend(loc=4) # Simple caption inside picture. 

### Plotting density evolution
ax2 = plt.subplot(plotB)
# plt.xlim((0,unitx*timePost.max()))
plt.xlim((xmin,xmax))
plt.ylim((ymin,ymax))

if(QuantityOfInterest == "density"):
  if(UnitSystem=="AU"):
    plt.ylabel('Density exc. el. (a.u.)')
  elif(UnitSystem=="SI"):
    plt.ylabel(r'Density exc. el. (m$^{-3}$)')
  plot_density, = ax2.plot(timeNe, observable_proj, 'k-', label=r'$n_e^{proj}(t)$')
  #ax2.legend(loc=3)
  if(PlottingMode=="keldysh"): #TODO: it is removed because too heavy. However, such information will be very interesting. 
    ax22 = ax2.twinx() #let's plot instant w_PI (t)
    #plot_wpi_av, = ax22.plot(timeNe[1:], wPI_TDDFT_average, 'b-', label=r'$<w_{PI}^{proj}>(t)$') #NOTHING PHYSICAL? 
    plot_wpi_approx, = ax22.plot(timeNe, wPI_TDDFT_approx, 'b-', label=r'$\widetilde{w_{PI}}=n_{exc}/2\tau$')
    if("STH" in PulseType or "StH" in PulseType):
      plot_wpi_Gruzdev, = ax22.plot(timeNe, PeakWpi_Gruzdev*np.ones(np.shape(timeNe)), 'b--', label=r'$w_{PI}$ [Gruzdev]') #NOTE: shows the peak value of w_pi 
  # if(EnableLogScale==1): 
  #  plt.plot(timeNe[1:], wPI_TDDFT_env.real, 'b--', label=r'$w_{PI}^{proj}(t)/cos(omega t)$')
    #plt.legend(loc=3)
    plt.ylabel('$w_{PI}(t)$')
    ax22.yaxis.label.set_color(plot_wpi_approx.get_color()) #colorizes the label
    ax22.spines["right"].set_edgecolor(plot_wpi_approx.get_color()) #colorizes the axis
    ax22.tick_params(axis='y', colors=plot_wpi_approx.get_color()) #colorizes the tics and numbers
  plt.xlim((xmin,xmax))
#  if(EnableLogScale==1): 
#   plt.ylim((ymin_wpi,ymax_wpi))
  #ax22.set_yscale('log')
#  plt.plot(unitx*timePost, observable_post, 'k--', label=r'TD-DFT: \frac{1}{V_{tot}} $\int_{\Omega} \left( n_{td}(r,t) - n_{gs}(r) \right)^{+} d^{3}r$')
#  plt.plot(timeKeldysh*unitx, N_Keldysh_SI, 'r-', label=r'[Keldysh (1964)]')
  ## when E=0, Keldysh is not defined. Then we have to restrict plotting where field is non zero.
  ## TODO: BUG: timeKeldysh and N_Gruzdev_SI have a different dimension! 
  
  ## If we plot anything else than STH, we can compute predictions of density using the Keldysh model. 
  if(PulseType == "SinSquared" or PulseType == "gaussian"): #then we plot performed predictions. 
    # Checking sizes of the data to plot. 
    if(len(timeKeldysh) != len(N_Gruzdev_SI)):
      print Header+"** Warning: tmin and tmax may be different in libKeldysh and in this file. "
      print Header+"Length(timeKeldysh) = "+str(len(timeKeldysh))
      print Header+"Length(N_Gruzdev_SI) = "+str(len(N_Gruzdev_SI))
    else: 
      timeKeldysh  = instants
      timeKeldyshB = instants
  
    if(PlottingMode=="doublepulse"): #NOTE: in ax2, plot only densities !
      ax2.plot(timeKeldysh*unitx,  N_Gruzdev_SI[:timeKeldysh.size] , 'b--', label=r'$n_{exc}[E_A(t)]$ [Gruzdev (2014)]')
      ax2.plot(timeKeldyshB*unitx, N_GruzdevB_SI[:timeKeldyshB.size], 'g--', label=r'$n_{exc}[E_B(t)]$ [Gruzdev (2014)]')
      ax2.plot(timeKeldysh*unitx,  N_Gruzdev_SI[:timeKeldysh.size]+N_GruzdevB_SI[:timeKeldysh.size], 'r--', label=r'$n_{exc}[E_A(t)] + n_{exc}[E_B(t)]$')
      #ax2.plot(timeKeldysh*unitx,  N_GruzdevAB_SI, 'r-', label=r'$n_{exc} [E_A(t)+E_B(t)]$') #BUG: UNCLEAR CASE
    else: #single pulse case
      ax2.plot(timeKeldysh*unitx,  N_Gruzdev_SI , 'k--', label=r'Exc. $n_{exc}(t)$ [Gruzdev (2014)]')
  #elif(PulseType=="STH"): #if STH, we only plot w_PI. 
    #ax2.plot(time*unitx,  PeakWpi_Gruzdev*np.ones(np.shape(time)), 'k--', label=r'Exc. $w_{PI}$ at max. [Gruzdev (2014)]')
    
elif (QuantityOfInterest == "number"):
  plt.ylabel('Number of exc. electrons')
  plt.plot(timeNe, Nproj_AU, 'k-', label=r'$N_e^{proj}(t)$')

#  plt.plot(unitx*timePost, Npost_AU, 'b-', label=r'$\int_{\Omega} \left( n_{td}(r,t) - n_{gs}(r) \right)^{+} d^{3}r$')
#plt.legend(loc=1)
# plt.legend(bbox_to_anchor=(CaptionX, CaptionY), loc=2, borderaxespad=0.)
if("STH" in PulseType or "StH" in PulseType):
    # Building the shared legend for subfig1 (if STH pulses)
    plot2   = [plot_density, plot_wpi_approx, plot_wpi_Gruzdev] #plot_wpi_av
    labels2 = [l.get_label() for l in plot2]
    ax2.legend(plot2, labels2, loc=2)

## Plotting vertical bars for comparison A(t) and Ne(t)
if(PlotGrid==1):
  for xc in xposition:
    plt.axvline(x=xc, color='#b3b3b3', linestyle=':')
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
plot_energy_electron, = ax3.plot(UnitTime * TimeAUtoSI * timeEnergies, Energy_Hartree_to_eV(ElectronEnergy), 'k-', label=r'$\xi_{el} (tot)$')
plot_energy_classical_potential, = ax3.plot(timeNe, (Egap_SI/e) * observable_proj * SuperCellVolume_SI, 'b-', label=r"$E_{gap} \times N_{exc}$ (eV)") #classical formula for electron gap transition 
plot_energy_classical_potential, = ax3.plot(timeNe, 2. * hbar * omega / e * observable_proj * SuperCellVolume_SI, 'b--', label=r"$2 \hbar \omega N_{exc}$ (eV)") #classical formula number 2 for 2-photon absorption #NOTE: wrong if MoveIons=Yes !
#plt.plot(UnitTime * TimeAUtoSI * timeEnergies, IonEnergy, 'b--', label=r"$\xi_{ion}$") #TODO: make automatic function of MoveIon value. 

# plt.plot(UnitTime * TimeAUtoSI * timeEnergies, Hartree, 'r-', label=r"$\xi_{Hartree}$")
plt.xlabel('Time (fs)')
plt.ylabel('Energy (eV)')
plt.legend(bbox_to_anchor=(CaptionX, CaptionY), loc=2, borderaxespad=0.)
plt.legend(loc=2)
## Plotting vertical bars for comparison A(t) and Ne(t)
if(PlotGrid==1):
  for xc in xposition:
    plt.axvline(x=xc, color='#b3b3b3', linestyle=':')

## ============ PLOTTING FOURIER SPECTRA =========
ax4 = plt.subplot(plotD)
# plt.plot(HHG_AU_order[:maxHHG], np.abs(HHG_AU[:maxHHG])**2, '-', label=r'HHG$(\omega)$')
plt.plot(HHG_AU_order[:maxHHG], np.power(HHG_AU_order[:maxHHG],2.)*np.abs(HHG_AU[:maxHHG])**2, '-', label=r'HHG$(\omega)$')
# ax4.set_xscale("log", nonposx='clip')
ax4.set_yscale("log", nonposx="clip")
plt.xlabel("Harmonic order")
plt.ylabel("Harmonic ampl. (arb. u.)")

## ======== Applying logscales...

if(EnableLogScale == 1):
  ax2.set_yscale("log", nonposy='clip')
  ax22.set_yscale("log", nonposy='clip')
  #ax2.legend(loc=4)
#  ax3.set_yscale("log", nonposy='clip')
#  ax3.legend(loc=4)
#else:
#  ax3.legend(loc=1)
  #ax2.legend(bbox_to_anchor=(CaptionX, CaptionY), loc=1, borderaxespad=0.)

## ==== Writing final picture =====

filename='KeldyshAnalysis-'+PlottingMode
if(EnableLogScale==1): 
  filename=filename+"-LogScale"
plt.savefig(filename+'.eps', bbox_inches='tight')
plt.savefig(filename+'.png', bbox_inches='tight')
if(ShowFinalPlot==1):
  plt.show()

print "** Info: Script ended sucessfully."
print "** Goodbye."
