#!/usr/bin/env python
#-*- coding: utf-8 -*-

# @module plotFinalQuantities
# Collects the generated octopus simulations from subfolders
# Builds a table with laser parameters, final densities, and final energies after simulation. 

## IMPORT LIBRARIES
import numpy as np
from numpy import genfromtxt, loadtxt, chararray
# from pylab import *
import sys
import decimal
import glob

# import calculateVolume as Sample
# import libKeldysh as Keldysh
from libAtomicUnits import *
from libPlottingSettings import *
from libString import *
from libPlotting import *

# List all output_iter files
for name in glob.glob('*/*/output_iter/'): 
  #TODO: this is not general enough, directory level does not help
  print name

  # Import the quantities from the subfolder
  geometry, Ne, laser, TotalCurrent, Energy = LoadRunOutputs(name)
  # Collect the laser parameters
  filename = name+"../inp"
  wavelength_inp = LookUpStringInFile(filename, "lambdaA =")
  wavelength_SI = float(wavelength_inp[10:17])
  PulseDuration_inp = LookUpStringInFile(filename, "tauAsi = ")
  PulseDuration_SI = float(PulseDuration_inp[8:15])
  DoublePulseDelay_inp = LookUpStringInFile(filename, "pulseDelaySI = ")
  DoublePulseDelay_SI = float(DoublePulseDelay_inp[15:22])
  
  print "Wavelength (nm): "+str(wavelength_SI)
  print "Pulse duration (s): "+str(PulseDuration_SI)
  print "Double pulse delay (s): "+str(DoublePulseDelay_SI)
  
  # Catch the peak field
  # /!\ TODO: duplicated code
  time = np.asfarray(laser[:,1]) #column 2
  A1x = np.asfarray(laser[:,2]) #column number 3
  A2x = np.asfarray(laser[:,5]) #column number 6
  timestepTable = np.diff(time) #size N-1
  TimeStep = timestepTable[0] #they are all equal, just select the first
  Efield = - (np.diff(A1x)/np.diff(time)) #correct value is *1/c 
  Efield2= - (np.diff(A2x)/np.diff(time))
  Efield_amplitude_SI=np.max(abs(Field_AU_to_SI(Efield / c_AU)))
  
  print "Electric field amplitude = "+str(Efield_amplitude_SI / 1E9)+" V/nm."
  
  del time, A1x, A2x, timestepTable, Efield, Efield2
  # Catch the final excited electron density. 
  
  # Append to the table

# Prepare the summary plots
## Plot as function of wavelength for one given peak intensity
## Plot as function of peak intensity for one given wavelength.
## Plot as function of double pulse delay


