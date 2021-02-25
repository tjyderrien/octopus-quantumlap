#!/usr/bin/env python2
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2020 T. J.-Y. Derrien
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

## @package libPulseGenerator
# This module defines various types of laser pulses and calls the spp-extended-theory/Keldysh/libKeldysh package to construct a pulse envelope for octopus. 
# A routine provides the pulse envelope in an octopus:td_from_file-compatible format. 

# IMPORT LIBRARIES
import numpy as np
from numpy import genfromtxt, loadtxt, chararray
from scipy.optimize import fsolve, root
import cmath
import sys
import matplotlib as mp
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
from matplotlib import rc
# from pylab import *
from scipy.constants import c, epsilon_0, mu_0, pi, e, m_e, h
from matplotlib.legend_handler import HandlerLine2D

# Settings for matplotlib
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':'16'})
## for Palatino and other serif fonts use:
rc('font',**{'family':'serif','serif':['Palatino'], 'size':'16'})
rc('text', usetex=True)
mp.rcParams['legend.numpoints'] = 1

## laser pulsation (SI system)
def omega(wavelength):#{{{
  return 2.0*pi*c/wavelength
#}}}

## normalization for gaussian width
def sigma(value): #{{ 
  return value/(2e0*cmath.sqrt(2e0*cmath.log(2e0)))
#}}

## 
def NonChirpedGaussianPulse(t, amplitude, wavelength, phase, duration, polarization, delay): #{{
  # Input: 1:R, 2:R, 3:R, 4:R, 5:R, 6:string, 7:R)
  # Output: R
  # defines the general temporal form of a non-chirped laser pulse
  #ElectricField = ((t-delay))**2e0 # 
  #TODO: generalize to non-normal incidence
  #TODO: include the polarization direction
  ElectricField = amplitude * cmath.cos(omega(wavelength)*t + phase) * cmath.exp(-0.5e0 * ((t - delay)/sigma(duration))**2e0)
  return ElectricField
#}}

NonChirpedGaussianPulse = np.vectorize(NonChirpedGaussianPulse)

## From [t, Pulse] table generates the figure of Pulse(t) in FigName file along with Title. 
def PlotPulse(t, Pulse, FigName, Title): #{{
  # Prepare some figures of required pulse durations
  maxy=3.0; miny = -maxy
  #Pulse = Pulse1+Pulse2
  Pulse = Pulse.real
  
  plt.figure()
  plt.xlabel('t $(fs)$')
  plt.ylabel('$E$ (V$/$nm)')
  plt.ylim(miny,maxy)
  plt.plot(1e15*t, Pulse*1E-9, label=Title)
  #plt.plot(t, Pulse1, label='Pulse 1')
  #plt.plot(t, Pulse2, label='Pulse 2')
  #plt.title(Title)
  plt.legend(loc=2)
  plt.tight_layout()
  plt.savefig(str(FigName)+'.eps')
  #plt.show()
#}}

## Print unrolled parameter list from a Parameter Table. 
def printCurrentParameters(ParameterTable): #{{
  duration1, amplitude1, wavelength1, phase1, polarization1, delay1 = ParameterTable[0, :]
  duration2, amplitude2, wavelength2, phase2, polarization2, delay2 = ParameterTable[1, :]
  
  tabValue=20
  print(("=== PULSE 1 ==="+"\t\t\t"+"=== PULSE 2 ===").expandtabs(tabValue))
  print(("QUANTITY"+"\t"+"VALUE"+"\t"+"UNIT"+"\t"+"QUANTITY"+"\t"+"VALUE"+"\t"+"UNIT").expandtabs(tabValue))
  print(("Wavelength1"+"\t"+str(wavelength1)+"\t"+" m"+"\t"+"Wavelength2"+"\t"+str(wavelength2)+"\t"+" m.").expandtabs(tabValue))
  print(("Phase1       "+"\t"+str(float(phase1)*180./pi)+"\t"+"(deg)"+"\t"+"Phase2"+"\t"+str(float(phase2)*180./pi)+"\t"+"(deg)").expandtabs(tabValue))
  print(("Polarization1"+"\t"+str(polarization1)+"\t"+"\t"+"Polarization2"+"\t"+str(polarization2)).expandtabs(tabValue))
  print(("Duration1"+"\t"+str(duration1)+"\t"+" s."+"\t"+"Duration2"+"\t"+str(duration2)+"\t"+" s.").expandtabs(tabValue))
  print(("Amplitude1"+"\t"+str(amplitude1)+"\t"+" V/m"+"\t"+"Amplitude2"+"\t"+str(amplitude2)+"\t"+" V/m.").expandtabs(tabValue))
  print(("Delay1"+"\t"+str(delay1)+"\t"+" s"+"\t"+"Delay2"+"\t"+str(delay2)+"\t"+" s").expandtabs(tabValue))
  print("")
  
#}}

## Extracts from a 2x6 data array the parameters of two laser pulses
def ExtractParameters(ParameterTable):#{{
  duration1, amplitude1, wavelength1, phase1, polarization1, delay1 = ParameterTable[0, :]
  duration2, amplitude2, wavelength2, phase2, polarization2, delay2 = ParameterTable[1, :]

  duration1 = np.asfarray(duration1)
  amplitude1 = np.asfarray(amplitude1)
  wavelength1 = np.asfarray(wavelength1)
  phase1 = np.asfarray(phase1)
  #polarization1 = np.asfarray(polarization1)
  delay1 = np.asfarray(delay1)

  duration2 = np.asfarray(duration2)
  amplitude2 = np.asfarray(amplitude2)
  wavelength2 = np.asfarray(wavelength2)
  phase2 = np.asfarray(phase2)
  #polarization1 = np.asfarray(polarization1)
  delay2 = np.asfarray(delay2)
  return duration1, amplitude1, wavelength1, phase1, polarization1, delay1, duration2, amplitude2, wavelength2, phase2, polarization2, delay2
#}}
