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

import calculateVolume as Sample
# import libKeldysh as Keldysh
from libAtomicUnits import *
from libPlottingSettings import *
from libString import *
from libPlotting import *

# List all output_iter files
for FolderName in glob.glob('*/*/output_iter/'): 
  #TODO: this is not general enough, directory level does not help
  print FolderName

  # Import the quantities from the subfolder
  geometry, Ne, laser, TotalCurrent, Energy = LoadRunOutputs(FolderName)
  # Collect the laser parameters
  filename = FolderName+"../inp"
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
  Efield_AU =Efield/c_AU #fixing dimension #TODO: clean me that...
  Efield2_AU=Efield2/c_AU
  Efield_amplitude_SI=np.max(abs(Field_AU_to_SI(Efield_AU)))
  
  print "Electric field amplitude = "+str(Efield_amplitude_SI / 1E9)+" V/nm."
  
  # Capturing cell volume
  SuperCellVolume_AU = Sample.VolumeOfSupercell_Log(FolderName+"../")
  SuperCellVolume_SI = Volume_AU_to_SI(SuperCellVolume_AU)
  SampleCellVolume_AU = SuperCellVolume_AU - Sample.VolumeOfVacuum(FolderName+"../") #Volume of the sample (without vacuum region) 
  CellVolume_AU = Sample.VolumePerPoint_Log(FolderName+"../") #0.06385 #0.07820 #Bohr**3
  SampleCellVolume_SI = Volume_AU_to_SI(SampleCellVolume_AU)
  print ""
  print "Super-cell total volume = "+str(SuperCellVolume_AU)+" Bohr**3."
  print "Super-cell total volume = "+str(SuperCellVolume_SI)+" m^{3}." 
  print ""
  print "Sample volume = "+str(SampleCellVolume_AU)+" Bohr**3."
  print "Sample vacuum volume = "+str(SampleCellVolume_SI)+" m^{3}."
  print ""
  print "Cell volume = "+str(CellVolume_AU)+" Bohr**3."
  print "Cell volume = "+str(Volume_AU_to_SI(CellVolume_AU))+" m^{-3}."

  # Catch the final excited electron density. 
  timeNe = np.asfarray(Ne[:,1])
  Nproj_AU = np.asfarray(Ne[:,2]) #number of particles, not a density!
  density_proj_AU = Nproj_AU / SampleCellVolume_AU #This is now a density (a.u.)
  
  ##==== Prepare plotting of energies ====
  timeEnergies   = np.asfarray(Energy[:,1]) #atomic units (! 1 is 2nd column)
  ElectronEnergy = np.asfarray(Energy[:,5]) #Hartree
  TotalEnergy    = np.asfarray(Energy[:,2]) #Hartree
  IonEnergy      = np.asfarray(Energy[:,3]) #Hartree
  Hartree        = np.asfarray(Energy[:,7]) #Hartree
  
  # Append to the table
  
  # Cleaning the guys before next loop
  del time, A1x, A2x, timestepTable, Efield, Efield2, geometry, Ne, laser, TotalCurrent, Energy
  

# Prepare the summary plots
## Plot as function of wavelength for one given peak intensity
## Plot as function of peak intensity for one given wavelength.
## Plot as function of double pulse delay


