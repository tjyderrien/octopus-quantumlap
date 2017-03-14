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
  
  # Append to the table

# Prepare the summary plots
## Plot as function of wavelength for one given peak intensity
## Plot as function of peak intensity for one given wavelength.
## Plot as function of double pulse delay


