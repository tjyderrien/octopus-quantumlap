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

# @module plotFinalQuantities
# Collects the generated octopus simulations from subfolders
# Builds a table with laser parameters, final densities, and final energies after simulation. 

## IMPORT LIBRARIES
import numpy as np
from numpy import genfromtxt, loadtxt, chararray
#from scipy import interpolate
import matplotlib.pyplot as pplt
# from pylab import *
import sys
import decimal

import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.legend_handler import HandlerLine2D
import mpl_toolkits.mplot3d as mplt 
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import interp2d

# import libKeldysh as Keldysh
#from libAtomicUnits import *
#from libPlottingSettings import *
from libPlotting import *
from libKeldysh import *
from libDatabase import *
import libExcitationDatabase

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':'14'})
## for Palatino and other serif fonts use:
#rc('font', **{'family':'serif', 'serif':['Palatino'], 'size':'16'})
rc('text', usetex=True) #True: Does not work on Draco. Just put False, then, and it will be smooth. 
mp.rcParams['legend.numpoints'] = 1

Header="[plotFinalQuantities] "

ShowFinalPlot              = 1

EnableDoublePulseAnalysis  = 2 #1: yes, 2: auto on command line
EnableInducedFieldAnalysis = 2 #1: yes, 2: auto on command line

Density=5E28 #specific to Si

if(EnableDoublePulseAnalysis==2 or EnableInducedFieldAnalysis==2): #automatic mode
  try: 
    EnableDoublePulseAnalysis = int(float(sys.argv[1]))
    EnableInducedFieldAnalysis = int(float(sys.argv[2]))
    print Header+"Input: Double pulse mode?"+str(EnableDoublePulseAnalysis)
    print Header+"Input: Induced fields   ?"+str(EnableInducedFieldAnalysis)
#    ShowFinalPlot = 0 #avoids problems with serial plotting 
  except: 
    print Header+"Usage: plotFinalQuantities.py DoublePulse? InducedFields? "
    print Header+"Ex:    plotFinalQuantities.py 1 1"
    exit()
  
# Capture the data
# For histocal reasons, two locations are possilbe.
path1='*/*/output_iter/' #valid for the runs computed from a same k-grid folder. 
path2='*/*/td.general/'

PlottingTable = libExcitationDatabase.GetExcitationDatabase(path2) #attempts to write a new DB. 
#if(len(PlottingTable) > 1): 
  #print Header+"** Info: capturing data from old PATH1." #BUG: destroys the previous dataase! 
#else:
  #PlottingTable = GetExcitationDatabase(path2)
  #print Header+"** Info: capturing data from novel path PATH2."

print PlottingTable

PlottingTable = FilterDatabase(PlottingTable, 0, IndexTable['FieldInduced']) #removes all induced fields data

runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingTable)

print Header+"** Preparing the summary plots..."


#### =============== PLOTTING STARTS HERE ==================
# The filtering will be done on the fly. Then we are sure of 
# what we plot in each figure. 

### Defining plotting structure

# Observables:        {N_exc^max, N_exc^final, irreversibility, Energy_max, Energy_final}

#### Plotting: parameters for controlling the canvas
if(EnableDoublePulseAnalysis==0):
  # 2 plots
  plotA=211; plotB=212; #plotC=313
elif (wavelength_numbers > 1):
  plotA=311; plotB=312; plotC=313;
else: #if(EnableDoublePulseAnalysis==1):
# 3 plots
# plotA=311; plotB=312; plotC=313;
  # 4 plots 
  plotA=221; plotB=222; plotC=223; plotD=224;

# Size of main canvas
SizeX = 16.
SizeY = 10.

#pplt.figure(figsize=(SizeX,SizeY))

## We change strategy: now we plot all the available data based on the number of parameters we have detected. 
# The input parameters will not define what we have to plot anymore. We prepare all the plots, which will simplify the understanding of this file. 
def plotEffectOfIntensity(PlottingTable, pulseDuration=10e-15, ShowFinalPlot=0):
  PlotVirtualExcitation = False
  print "#====== Plot: SINGLE PULSE, VARIOUS WAVELENGTHS, SAME PULSE DURATION, SAME POLARIZATIONS: final density as function of field ======"
  #Filter: double pulses, same wavelengths, same polarization
  PlottingSinglePulses = PlottingTable #get the full database
  PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
  PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0. 
  #PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['CEP1']) #removes fields where CEP != 0. #NOTE: Can be done for Si-LDA. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['FreezeHXC']) #removes fields where FreezeHXC != 0. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['FieldInduced']) #removes fields where FreezeHXC != 0. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 1., IndexTable['TheoryLevel']) #removes runs where TheoryLevel != dft. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 2, IndexTable['PulseType']) #keeps only pulses of type STH (2)
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, pulseDuration, IndexTable['pulses1']) #keeps only pulses of type STH (2)

  ## prints the current state of selection
  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq, wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = CountTheNumberOfLaserParameters(PlottingSinglePulses)

  # Simple pulse, single color
  pplt.figure()
  #if(len(wavelength_uniq) == 1): #{{{ #single color simulations
  # Complexity level: 1
  
  ax1 = pplt.subplot(plotA)
  
  pplt.title('Density as function of peak intensity')
  # Number of sweeped parameters: 3 (Efield, Delay, CrossPolarization)
  #pplt.scatter(FieldToIntensity(Efield1), Ne_final, c=Delay, s=(30.+CrossPolarization*50.), label=r'$N_{exc}$[Field]')
  for wavelength_value in wavelength_uniq: #{{{
    # select the corresponding wavelength in results
    PlottingSinglePulses_wavelength = FilterDatabase(PlottingSinglePulses, wavelength_value, IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts.

    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingSinglePulses_wavelength)

    pplt.plot(FieldToIntensity(Efield1), Ne_final, 'o', label=r'$N_{exc}$, $\lambda=$'+str(round(wavelength_value*1E9))+' nm')
  #}}}
  ax1.legend(loc=2, borderaxespad=0.) #bbox_to_anchor=(1.05, 1), 
  if(not PlotVirtualExcitation): 
    pplt.xlabel(r'Peak intensity (W/m$^2$)')
  pplt.ylabel(r'Final exc. density (m$^{-3}$)')
#  pplt.ylim((1E16,4.*5e28))
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  
  ax12 = ax1.twinx()
  for wavelength_value in wavelength_uniq: #{{{
    # select the corresponding wavelength in results
    PlottingSinglePulses_wavelength = FilterDatabase(PlottingSinglePulses, 800e-9, IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts.

    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingSinglePulses_wavelength)

    pplt.plot(FieldToIntensity(Efield1), Ne_final/pulses1, 'o') #, label=r'$N_{exc}/\tau$, $\lambda=$'+str(round(wavelength_value*1E9))+' nm')
  #}}}
  ax12.set_yscale('log')
  pplt.ylabel('$w_{PI}$ (m$^{-3}$.s$^{-1}$)')
  if(PlotVirtualExcitation):
    ax2 = pplt.subplot(plotB)
    pplt.title('Virtual exc. factor')
    pplt.plot(FieldToIntensity(Efield1), irreversibility, 'ob', label='Virtual exc.')
    pplt.xlabel(r'Intensity (W/m$^2$)')
    pplt.ylabel('Virtual exc. factor')
    pplt.ylim((0,1))
    ax2.set_xscale('log')
    ax22 = ax2.twinx()
    for wavelength_value in wavelength_uniq:
      # select the corresponding wavelength in results
      PlottingSinglePulses_wavelength = FilterDatabase(PlottingSinglePulses, wavelength_value, IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts.

      runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingSinglePulses_wavelength)
      pplt.plot(FieldToIntensity(Efield1), gammaKeldysh(2.56*e, 0.228, Efield1, wavelengths1), 'or', label=r'Adiabaticity, $\lambda=$'+str(round(wavelength_value*1E9))+' nm')
    pplt.ylabel('Adiabadicity factor')
    #pplt.ylim((0,1))
    ax22.set_yscale('log')
      
    # ax2.set_yscale("log")
    ax2.legend(loc=2, borderaxespad=0.)
    # ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax22.legend(loc=1, borderaxespad=0.)
  ax12.legend(loc=4, borderaxespad=0.)
  
  # Saving the plot
  ExportName = "EffectOfIntensity-"+str(pulseDuration*1E15)+"fs"
  pplt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.eps')
  pplt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.png')
  print Header+"** Saved: Metaplot-SinglePulse-SingleColor-"+ExportName+".eps|png."
  if(ShowFinalPlot==1):
    pplt.show()

  del ax1, ax12, PlottingSinglePulses, wavelength_uniq, CEP1_uniq
  if(PlotVirtualExcitation): 
    del ax2, ax22
#}}}

## We change strategy: now we plot all the available data based on the number of parameters we have detected. The input parameters will not define what we have to plot anymore. We prepare all the plots, which will simplify the understanding of this file. 
def plotEffectOfIntensity_PulseDuration(PlottingTable, pulseWavelength=800e-9):
  PlotVirtualExcitation = False
  print "#====== Plot: SINGLE PULSE, SAME WAVELENGTHS, SAME POLARIZATIONS: final density as function of field and pulse duration ======"
  #Filter: double pulses, same wavelengths, same polarization
  PlottingSinglePulses = PlottingTable #get the full database
  PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
  PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0. 
  #PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['CEP1']) #removes fields where CEP != 0. #NOTE: Can be done for Si-LDA. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['FreezeHXC']) #removes fields where FreezeHXC != 0. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['FieldInduced']) #removes fields where FreezeHXC != 0. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 1., IndexTable['TheoryLevel']) #removes runs where TheoryLevel != dft. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 2, IndexTable['PulseType']) #keeps only pulses of type STH (2)
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, pulseWavelength, IndexTable['wavelengths1']) #keeps only pulses of type STH (2)

  ## prints the current state of selection
  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = CountTheNumberOfLaserParameters(PlottingSinglePulses)

  # Simple pulse, single color
  pplt.figure()
  #if(len(wavelength_uniq) == 1): #{{{ #single color simulations
  # Complexity level: 1
  
  ax1 = pplt.subplot(plotA)
  
  pplt.title('Density as function of peak intensity')
  # Number of sweeped parameters: 3 (Efield, Delay, CrossPolarization)
  #pplt.scatter(FieldToIntensity(Efield1), Ne_final, c=Delay, s=(30.+CrossPolarization*50.), label=r'$N_{exc}$[Field]')
  for pulse_value in pulses_uniq: #{{{
    # select the corresponding wavelength in results
    PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses, pulse_value, IndexTable['pulses1'])
    PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses_pulses, pulse_value, IndexTable['pulses2']) #keeps only pulses of the same duration.

    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingSinglePulses_pulses)

    pplt.plot(FieldToIntensity(Efield1), Ne_final, 'o', label=r'$N_{exc}$, $\tau_1=$'+str(round(pulse_value*1E15))+r'fs, $\tau_2=$'+str(round(pulse_value*1E15))+'fs')
  #}}}
  ax1.legend(loc=2, borderaxespad=0.) #bbox_to_anchor=(1.05, 1), 
  #pplt.xlabel(r'Peak intensity (W/m$^2$)')
  pplt.ylabel(r'Final exc. density (m$^{-3}$)')
#  pplt.ylim((1E16,4.*5e28))
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  
  ax12 = ax1.twinx()
  for pulse_value in pulses_uniq: #{{{
    # select the corresponding wavelength in results
    PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses, pulse_value, IndexTable['pulses1'])
    PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses_pulses, pulse_value, IndexTable['pulses2']) #keeps only pulses of the same duration.

    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingSinglePulses_pulses)

    pplt.plot(FieldToIntensity(Efield1), Ne_final/pulses1, 'o') #, label=r'$N_{exc}/\tau$, $\lambda=$'+str(round(wavelength_value*1E9))+' nm')
  #}}}
  ax12.set_yscale('log')
  pplt.ylabel('$w_{PI}$ (m$^{-3}$.s$^{-1}$)')
  if(PlotVirtualExcitation):
    ax2 = pplt.subplot(plotB)
    pplt.title('Virtual exc. factor')
    pplt.plot(FieldToIntensity(Efield1), irreversibility, 'ob', label='Virtual exc.')
    pplt.xlabel(r'Intensity (W/m$^2$)')
    pplt.ylabel('Virtual exc. factor')
    pplt.ylim((0,1))
    ax2.set_xscale('log')
    ax22 = ax2.twinx()
    for pulse_value in pulses_uniq:
      # select the corresponding wavelength in results
      PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses, pulse_value, IndexTable['pulses1']) #keeps only 800 nm wavelenghts.

      runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingSinglePulses_pulses)
      pplt.plot(FieldToIntensity(Efield1), gammaKeldysh(2.56*e, 0.228, Efield1, wavelengths1), 'or', label=r'Adiabaticity, $\tau=$'+str(round(pulse_value*1E15))+' fs')
    pplt.ylabel('Adiabadicity factor')
    #pplt.ylim((0,1))
    ax22.set_yscale('log')
      
    # ax2.set_yscale("log")
    ax2.legend(loc=2, borderaxespad=0.)
    # ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax22.legend(loc=1, borderaxespad=0.)
  ax12.legend(loc=4, borderaxespad=0.)
  
  # Saving the plot
  ExportName = "EffectOfPulseIntensity-Duration-"+str(pulseWavelength*1E9)+"nm"
  pplt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.eps')
  pplt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.png')
  print Header+"** Saved: 'Metaplot-SinglePulse-SingleColor-'"+ExportName+".png."
  #if(ShowFinalPlot==1):
    #pplt.show()

  del ax1, ax12, PlottingSinglePulses, wavelength_uniq, CEP1_uniq
  if(PlotVirtualExcitation): 
    del ax2, ax22
#}}}

def plotEffectOfCEP(PlottingTable):
  print "#====== Plot: SINGLE PULSE, EFFECT OF THE CEP, SAME WAVELENGTHS, SAME POLARIZATIONS: final density as function of field ======"
  #Filter: double pulses, same wavelengths, same polarization

  # Simple pulse, single color
  PlottingSinglePulses = PlottingTable #get the full database
  PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
  PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 800e-9, IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts. #TODO: this field should be making a loop in available wavelengths given in wavelengths1_uniq.
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0, IndexTable['FreezeHXC']) #keeps only non-frozen XC. 
  #PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['CEP1']) #removes fields where CEP != 0.

  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = CountTheNumberOfLaserParameters(PlottingSinglePulses)
  
  # Complexity level: 2

  pplt.figure(figsize=(SizeX,SizeY))
  ax1 = pplt.subplot(111) #111 or plotA if using plotB
  pplt.title('Final density as function of peak intensity')
  # Number of sweeped parameters: 3 (Efield, Delay, CrossPolarization)
  #pplt.scatter(FieldToIntensity(Efield1), Ne_final, c=Delay, s=(30.+CrossPolarization*50.), label=r'$N_{exc}$[Field]')
  for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
    PlottingSinglePulses_CEP = FilterDatabase(PlottingSinglePulses, CEP1_value, IndexTable['CEP1']) #select the local CEP
    
    ## prints the current state of selection and stores the values of sweeped parameters
    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingSinglePulses_CEP)
    
    pplt.plot(FieldToIntensity(Efield1), Ne_final, 'o', label=r'$N_{exc}^{final}$, $\phi_{CE}$='+str(round(CEP1_value*180./pi)))
    
    ax1.legend(loc=4, borderaxespad=0.) #bbox_to_anchor=(1.05, 1), 
    #pplt.xlabel(r'Peak intensity (W/m$^2$)')
    pplt.ylabel(r'Final density (m$^{-3}$)')
    #  pplt.ylim((1E16,4.*5e28))
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    print Efield1, CEP1_value
    
  #Adding axis of gammaKeldysh, but for all CEPs !!
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingSinglePulses)
  
  ax12 = ax1.twinx()
  
  #pplt.plot(FieldToIntensity(Efield1), Ne_final/(pulses1), 'o') #, label=r'$N_{exc}/\tau$') #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
  
  ax12.set_yscale('log')
  ax12.set_ylabel(r'$w_{PI}$ (m$^{-3}$.s$^{-1}$)')

  #ax2 = pplt.subplot(plotB)
  #pplt.title('Virtual excitation factor')
  
  #pplt.plot(FieldToIntensity(Efield1), irreversibility, 'ow', label='Virtual exc.')
  
  #pplt.xlabel(r'Peak intensity (W/m$^2$)')
  #pplt.ylabel('Virtual exc. factor')
  #pplt.ylim((0,1))
  #ax2.set_xscale('log')
  ##adding the Keldysh adiabatic factor
  #ax22 = ax2.twinx()
  
  #pplt.plot(FieldToIntensity(Efield1), gammaKeldysh(2.56*e, 0.228, Efield1, wavelengths1), 'or', label='Adiab. fact.') #BUG: 2.56 must not be specified by hand here! 

  #pplt.ylabel('Adiabadicity factor')
  ##pplt.ylim((0,1))
  #ax22.set_yscale('log')

  #ax2.set_yscale("log")
  #ax2.legend(loc=2, borderaxespad=0.)
  ax1.legend(bbox_to_anchor=(1.15, 1), loc=2, borderaxespad=0.)
  #ax22.legend(loc=1, borderaxespad=0.)

  # Saving the plot
  ExportName = "EffectOfCEP"
  pplt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.eps')
  pplt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.png')
  #if(ShowFinalPlot==1):
    #pplt.show()

  del ax1, PlottingSinglePulses, wavelength_uniq, wavelengths1
  #}}}
 
def plotEffectOfFrozenXC(PlottingTable):
  print "#====== Plot: SINGLE PULSE, EFFECT OF THE FROZEN HXC, SAME WAVELENGTHS, SAME POLARIZATIONS: final density as function of field ======"
  #Filter: double pulses, same wavelengths, same polarization

  # Simple pulse, single color
  PlottingSinglePulses = PlottingTable #get the full database
  PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
  PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 800e-9, IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts. #TODO: this field should be making a loop in available wavelengths given in wavelengths1_uniq. 
  #PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0, IndexTable['CEP1']) #keeps only CEP=0. 

  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = CountTheNumberOfLaserParameters(PlottingSinglePulses)

  #if(len(CEP1_uniq) == 1): #{{{ #single color simulations
  # Complexity level: 2

  pplt.figure()
  ax1 = pplt.subplot(plotA)
  pplt.title('Final density as function of peak intensity')
  # Number of sweeped parameters: 3 (Efield, Delay, CrossPolarization)
  #pplt.scatter(FieldToIntensity(Efield1), Ne_final, c=Delay, s=(30.+CrossPolarization*50.), label=r'$N_{exc}$[Field]')
  #for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
  #PlottingSinglePulses_CEP = FilterDatabase(PlottingSinglePulses, CEP1_value, IndexTable['CEP1']) #select the local CEP
  
  # Filtering the FreezeHXC states
  PlottingSinglePulses_CEP_NoFreezeHXC = FilterDatabase(PlottingSinglePulses, 0, IndexTable['FreezeHXC']) #select the non-frozen HXC
  PlottingSinglePulses_CEP_FreezeHXC = FilterDatabase(PlottingSinglePulses, 1, IndexTable['FreezeHXC']) #select the frozen HXC
  ## prints the current state of selection and stores the values of sweeped parameters
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingSinglePulses_CEP_NoFreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), Ne_final, 'ow', label=r'$N_{exc}^{final}$') #, CEP='+str(round(CEP1_value*180./pi)))
  
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingSinglePulses_CEP_FreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), Ne_final, 'ok', label=r'$N_{exc}^{final}$, FrozenHXC')#, CEP='+str(round(CEP1_value*180./pi)))
  
  ax1.legend(loc=4, borderaxespad=0.) #bbox_to_anchor=(1.05, 1), 
  #pplt.xlabel(r'Peak intensity (W/m$^2$)')
  pplt.ylabel(r'Final density (m$^{-3}$)')
  #  pplt.ylim((1E16,4.*5e28))
  ax1.set_xscale('log')
  ax1.set_yscale('log')

  #Adding axis of gammaKeldysh
  ax12 = ax1.twinx()
  
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingSinglePulses_CEP_NoFreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), Ne_final/pulses1, 'ow') #, label=r'$N_{exc}/\tau$') #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
  
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingSinglePulses_CEP_FreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), Ne_final/pulses1, 'ok') #, label=r'$N_{exc}/\tau$, FrozenHXC') #, I='+str(FieldToIntensity(Efield1_value))+r'
  
  ax12.set_yscale('log')
  ax12.set_ylabel(r'$w_{PI}$ (m$^{-3}$/s)')

  ax2 = pplt.subplot(plotB)
  pplt.title('Virtual excitation factor')
  
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingSinglePulses_CEP_NoFreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), irreversibility, 'ow', label='Virtual exc.')
  
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingSinglePulses_CEP_FreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), irreversibility, 'ok', label='Virtual exc., FreezeHXC')
  
  pplt.xlabel(r'Peak intensity (W/m$^2$)')
  pplt.ylabel('Virtual exc. factor')
  pplt.ylim((0,1))
  ax2.set_xscale('log')
  #adding the Keldysh adiabatic factor
  ax22 = ax2.twinx()
  
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingSinglePulses_CEP_NoFreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), gammaKeldysh(2.56*e, 0.228, Efield1, wavelengths1), 'or', label='Adiabaticity') #BUG: 2.56 must not be specified by hand here! 
  
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingSinglePulses_CEP_FreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), gammaKeldysh(2.56*e, 0.228, Efield1, wavelengths1), 'or') #, label='Adiabaticity') #BUG: 2.56 must not be specified by hand here! 
  
  pplt.ylabel('Adiabadicity factor')
  #pplt.ylim((0,1))
  ax22.set_yscale('log')

  #print Efield1, CEP1_value
  # ax2.set_yscale("log")
  #ax2.legend(loc=2, borderaxespad=0.)
  # ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
  #ax22.legend(loc=1, borderaxespad=0.)

  # Saving the plot
  pplt.savefig('Metaplot-SinglePulse-SingleColor-EffectOfFrozenHXC.eps')
  pplt.savefig('Metaplot-SinglePulse-SingleColor-EffectOfFrozenHXC.png')
  #if(ShowFinalPlot==1):
    #pplt.show()

  del ax1, PlottingSinglePulses, wavelength_uniq, wavelengths1
  #}}}

#del PlottingTable

## Plots the effect of wavelength (1D plot) by introducing 10 percents tolerance on the electric field
def plotEffectOfWavelength_Tolerance(PlottingTable):
  print "#====== Plot: SINGLE PULSE, EFFECT OF THE WAVELENGTH, SAME CEP, SAME POLARIZATIONS: final density as function of field ======"
  #Filter: double pulses, same wavelengths, same polarization

  # Simple pulse, single color
  PlottingSinglePulses = PlottingTable #get the full database
  PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
  PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['CrossPolarization']) #removes fields where CrossPolarization is different from 0. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['CEP1']) #keeps only CEP = 0. 

  ## TODO: check if CEP could count at smaller wavelength. 

  print Header+"Available data: "
  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = CountTheNumberOfLaserParameters(PlottingSinglePulses)
  #if(len(CEP1_uniq) == 1): #{{{ #single color simulations
  # Complexity level: 2

  pplt.figure(figsize=(SizeX,SizeY)) #new figure
  ax1 = pplt.subplot(111)
  pplt.title('Final density as function of laser wavelength')
  print Header+"Total number of datasets: "+str(len(set(Efield1_uniq)))
  for Efield1_value in Efield1_uniq: #for each value of Efield1, we make a new search, and plot one new line in the chart
    #PlottingSinglePulses_Efield = FilterDatabase(PlottingSinglePulses, Efield1_value, IndexTable['Efield1']) #select the local Efield
    # TODO: how to merge data with values close to each others? Maybe introduce tolerance of 10% ?
    PlottingSinglePulses_Efield = FilterDatabaseLowerThan(PlottingSinglePulses, Efield1_value*1.1, IndexTable['Efield1'])
    PlottingSinglePulses_Efield = FilterDatabaseGreaterThan(PlottingSinglePulses_Efield, Efield1_value*0.9, IndexTable['Efield1'])
    print Header+"Dataset for the field: "+str(Efield1_value*1e-9)+' V/nm'
    #print PlottingSinglePulses_Efield
    # TODO: as fields are all slighly different, we should project on the neareast neighboring value, up to 1E8 V/m. 
    ## prints the current state of selection and stores the values of sweeped parameters
    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingSinglePulses_Efield)
    
    #print Efield1_value
    # TODO: PUT here a select case which attributes a new color to a new field interval
    if(0.99e8 < Efield1_value < 1.01e8):
      pplt.plot(1e9*wavelengths1, Ne_final, 'or', label=r'$N_{exc}^{final}$') #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
    if(0.99e9 < Efield1_value < 1.01e9):
      pplt.plot(1e9*wavelengths1, Ne_final, 'ob', label=r'$N_{exc}^{final}$') #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
    if(4.99e9 < Efield1_value < 5.01e9):
      pplt.plot(1e9*wavelengths1, Ne_final, 'og', label=r'$N_{exc}^{final}$') #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
    if(8.99e9 < Efield1_value < 9.01e9):
      pplt.plot(1e9*wavelengths1, Ne_final, 'ok', label=r'$N_{exc}^{final}$') #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
    del PlottingSinglePulses_Efield
      
  #ax1.legend(loc=3, borderaxespad=0., bbox_to_anchor=(1.05, 0.05))
  pplt.xlabel(r'Laser wavelength (nm)')
  pplt.ylabel(r'Final density (m$^{-3}$)')
  pplt.xlim((100.,1500.))
  #ax1.set_xscale('log')
  ax1.set_yscale('log')
  ax12 = ax1.twinx()
  for Efield1_value in Efield1_uniq: #for each value of Efield1, we make a new search, and plot one new line in the chart
    #PlottingSinglePulses_Efield = FilterDatabase(PlottingSinglePulses, Efield1_value, IndexTable['Efield1']) #select the local Efield
    # TODO: how to merge data with values close to each others? Maybe introduce tolerance of 10% ?
    PlottingSinglePulses_Efield = FilterDatabaseLowerThan(PlottingSinglePulses, Efield1_value*1.1, IndexTable['Efield1'])
    PlottingSinglePulses_Efield = FilterDatabaseGreaterThan(PlottingSinglePulses_Efield, Efield1_value*0.9, IndexTable['Efield1'])
    print Header+"Dataset for the field: "+str(Efield1_value*1e-9)+' V/nm'
    #print PlottingSinglePulses_Efield
    # TODO: as fields are all slighly different, we should project on the neareast neighboring value, up to 1E8 V/m. 
    ## prints the current state of selection and stores the values of sweeped parameters
    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingSinglePulses_Efield)
    
    #print Efield1_value
    # TODO: PUT here a select case which attributes a new color to a new field interval
    if(0.99e8 < Efield1_value < 1.01e8):
      pplt.plot(1e9*wavelengths1, Ne_final/pulses1, 'or', label=r'$N_{exc}/\tau$') #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
    if(0.99e9 < Efield1_value < 1.01e9):
      pplt.plot(1e9*wavelengths1, Ne_final/pulses1, 'ob', label=r'$N_{exc}/\tau$') #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
    if(4.99e9 < Efield1_value < 5.01e9):
      pplt.plot(1e9*wavelengths1, Ne_final/pulses1, 'og', label=r'$N_{exc}/\tau$') #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
    if(8.99e9 < Efield1_value < 9.01e9):
      pplt.plot(1e9*wavelengths1, Ne_final/pulses1, 'ok', label=r'$N_{exc}/\tau$') #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
    del PlottingSinglePulses_Efield
  ax12.set_yscale('log')
  ax12.set_ylabel(r'$w_{PI}$ (m$^{-3}$/s)')

  #ax2 = plt.subplot(plotB)
  #plt.title('Irreversibility factor')
  #pplt.plot(FieldToIntensity(Efield1), irreversibility, 'o', label='Irreversibility')
  #pplt.xlabel(r'Intensity (W/m$^2$)')
  #pplt.ylabel('Irreversibility factor')
  #pplt.ylim((0,1))
  #ax2.set_xscale('log')
  #ax22 = ax2.twinx()
  #pplt.plot(FieldToIntensity(Efield1), gammaKeldysh(2.56*e, 0.228, Efield1, wavelengths1), 'or', label='Adiabaticity')
  #pplt.ylabel('Adiabadicity factor')
  ##pplt.ylim((0,1))
  #ax22.set_yscale('log')

  # ax2.set_yscale("log")
  #ax2.legend(loc=2, borderaxespad=0.)
  # ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
  #ax22.legend(loc=1, borderaxespad=0.)

  # Saving the plot
  pplt.savefig('Metaplot-SinglePulse-SingleColor-EffectOfWavelength_Tol10.eps')
  pplt.savefig('Metaplot-SinglePulse-SingleColor-EffectOfWavelength_Tol10.png')
  #if(ShowFinalPlot==1):
    #pplt.show()

  del ax1, PlottingSinglePulses
  #}}}

## Plots the N_exc(wavelength, field) (2D plot)
def plotEffectOfWavelength_2dMap(PlottingTable):
  print "#====== Plot: SINGLE PULSE, EFFECT OF THE WAVELENGTH, SAME CEP, SAME POLARIZATIONS: final density as function of field ======"
  #Filter: double pulses, same wavelengths, same polarization

  # Simple pulse, single color
  PlottingSinglePulses = PlottingTable #get the full database
  PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
  PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['CrossPolarization']) #removes fields where CrossPolarization is different from 0. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['CEP1']) #keeps only CEP = 0. 

  ## TODO: check if CEP could count at smaller wavelength. 

  print Header+"Available data: "
  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = CountTheNumberOfLaserParameters(PlottingSinglePulses)
  
  #plt.figure(figsize=(SizeX,SizeY)) #new figure
  DensityMax2D = np.empty((0,3)) #preparing the 2D plot
  #print DensityMax2D
  #ax1 = pplt.subplot(111)
  plt.title('Final density as function of laser wavelength')
  print Header+"Total number of datasets: "+str(len(set(Efield1_uniq)))
  
  ##Building up the stack to plot here, therefore we could plot anytype of plot. 
  for Efield1_value in Efield1_uniq: #for each value of Efield1, we make a new search, and plot one new line in the chart
    PlottingSinglePulses_Efield = FilterDatabase(PlottingSinglePulses, Efield1_value, IndexTable['Efield1']) #select the local Efield
    for wavelength_value in wavelength_uniq:
      PlottingSinglePulses_Efield_wavelength = FilterDatabase(PlottingSinglePulses_Efield, wavelength_value, IndexTable['wavelengths1']) #select the local Efield
      runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingSinglePulses_Efield_wavelength)
      
      if(len(Ne_final) > 0e0): #Filtering system is weird
        #DensityMax2D=np.vstack((DensityMax2D, [wavelengths1[0], Efield1[0], Ne_final[0]]))
        ToBeAdded=np.array([wavelengths1[0], Efield1[0], Ne_final[0]])
        #print "Stacking..."
        #print ToBeAdded
        DensityMax2D=np.vstack((DensityMax2D, ToBeAdded))
  
  # Stack is built. Export, show and plot on a heat map. 
  print "** Exporting the sub-plot2D: Ne_final(wavelength, Field)... Useful to visualize available data with GNUplot. "
  ExportToTxt(DensityMax2D, "Wavelength-Field-NeFinal.csv")
  plt.figure()
  #CS = plt.imshow(DensityMax2D, cmap=plt.cm.Blues) #WORKS, but data are too small (no interpolation)
  CS = plt.contourf(DensityMax2D, cmap=plt.cm.Blues) #WORKS, but there are degenerated data + not enough data for interpolation. 
  #CS = plt.contourf(1E9*np.transpose(DensityMax2D[:,0]),1E9*DensityMax2D[:,1],np.log10(DensityMax2D), cmap=plt.cm.RdBu_r) #WORKS
  plt.colorbar(CS)
  #plt.contourf(wavelength_uniq, Efield1_uniq, DensityMax2D) #, label=r'$N_{exc}^{final}$') #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
  del PlottingSinglePulses_Efield
      
  #ax1.legend(loc=3, borderaxespad=0., bbox_to_anchor=(1.05, 0.05))
  plt.xlabel(r'Laser wavelength (um)')
  plt.ylabel(r'Field amplitude (cell number)')
  #plt.ylabel(r'Final density (m$^{-3}$)')
  #plt.ylabel(r'Final density (m$^{-3}$)')
  #plt.xlim((100e-9,3100e-9))
  #ax1.set_xscale('log')
  #ax1.set_yscale('log')

  # Saving the plot
  pplt.savefig('Metaplot-SinglePulse-SingleColor-EffectOfWavelengthAndField.eps')
  pplt.savefig('Metaplot-SinglePulse-SingleColor-EffectOfWavelengthAndField.png')
  #if(ShowFinalPlot==1):
    #pplt.show()

  del PlottingSinglePulses, DensityMax2D
#}}}


## Orgazing the complicated life of double-pulse plotting. 
# Mode 1: fixed wavelengths, same pulse duration, then plot Ne_final(Efield1, Efield2)
def plotEffectOfDoublePulse_Efield1_Efield2_batchAllColors(PlottingTable, PulseDuration=20e-15): #{{{
  
  print "====== Plot: DOUBLE PULSE, SAME POLARIZATION, SAME PULSE DURATIONS ====="
  #Filters: double pulses, same wavelengths, same polarization
  PlottingDoublePulses = PlottingTable #get the full database
  PlottingDoublePulses = FilterDatabaseGreaterThan(PlottingDoublePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak < 1 V/m. #works
  PlottingDoublePulses = FilterDatabaseGreaterThan(PlottingDoublePulses, 1E0, IndexTable['Efield2']) #removes fields where Epeak2 < 1 V/m to select double pulses.
  PlottingDoublePulses = FilterDatabase(PlottingDoublePulses, float(0.), IndexTable['CrossPolarization']) #removes fields where Polarization1 != Polarization2
  PlottingDoublePulses = FilterDatabase(PlottingDoublePulses, PulseDuration, IndexTable['pulses2']) # Removing fields were pulseDuration1 != pulseDuration2
  
  #### 2D PLOTTING Density(Efield1, Efield2) for one given set of laser parameters. For each combination, we produce one plot like this. 
  
  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq, wavelength2_uniq, CEP2_uniq, Efield2_uniq, pulses2_uniq= CountTheNumberOfLaserParameters(PlottingDoublePulses)
  #plt.figure()
  print Header+"Total number of datasets: "+str(len(set(wavelength_uniq))*len(set(wavelength2_uniq)))
  for Color1 in wavelength_uniq: #selecting color1
    PlottingDoublePulses_SelectedColor1 = FilterDatabase(PlottingDoublePulses, Color1, IndexTable['wavelengths1'])
    for Color2 in wavelength2_uniq: #selecting color2
      # Selecting Color 2
      PlottingDoublePulses_SelectedColor2 = FilterDatabase(PlottingDoublePulses_SelectedColor1, Color2, IndexTable['wavelengths2'])
      
      # Constructing the density table
      DensityMax2D = np.empty((0,3)) #preparing the 2D plot
      runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingDoublePulses_SelectedColor2)
      
      print "Color1: "+str(Color1)+" m, Color2: "+str(Color2)+" m."
      if(len(Ne_final) >= 1):
        ToBeAdded=np.array([Efield1[0], Efield2[0], Ne_final[0]]) #BUG: removes some results, not clean enough!
        print ToBeAdded
        DensityMax2D=np.vstack((DensityMax2D, ToBeAdded))
        print "Cleaning..."
        del ToBeAdded
        # Stack is built. Export, show and plot on a heat map. 
        print "** Info: exporting the sub-plot2D: Ne_final(wavelength, Field)..."
        print "** Info: useful to visualize available data with GNUplot. "
        ExportToTxt(DensityMax2D, "Wavelength1-"+str(round(Color1*1E9))+"-Wavelength2-"+str(round(Color2*1E9))+"-NeFinal.csv")
      
      # If the size is big enough, let's make a plot.
      if(len(Ne_final) >= 2):
        plt.figure()
        plt.title(r'$N_{exc}(E_1,E_2)$, $\lambda_1=$'+str(Color1*1E9)+r'nm, $\lambda_2=$'+str(Color2*1E9)+'nm.')
        CS = plt.contourf(DensityMax2D, cmap=plt.cm.Blues)
        plt.colorbar(CS)
      
        filename='Metaplot-DoublePulse-'+str(round(Color1*1E9))+'nm-'+str(round(Color2*1E9))+'nm'
        plt.savefig(filename+'.eps')
        plt.savefig(filename+'.png')
        if(ShowFinalPlot==1):
          plt.show()
      
      del DensityMax2D 
#}}}
  


##### 2D PLOTTING Density(wavelength1, wavelength2) for one single field value? 
#Mode 2: fixed fields, same pulse durations, then plot Ne_final(Color1, Color2)
def plotEffectOfDoublePulse_Color1_Color2_SingleFieldStrength(PlottingTable, FieldStrength, PulseDuration=20e-15, Tolerance=0.01): #{{{
  print "====== Plot: DOUBLE PULSE, SAME POLARIZATION, SAME PULSE DURATIONS ====="
  #Filters: double pulses, same wavelengths, same polarization
  PlottingDoublePulses = PlottingTable #get the full database
  PlottingDoublePulses = FilterDatabaseGreaterThan(PlottingDoublePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak < 1 V/m. #works
  PlottingDoublePulses = FilterDatabaseGreaterThan(PlottingDoublePulses, 1E0, IndexTable['Efield2']) #removes fields where Epeak2 < 1 V/m to select double pulses.
  PlottingDoublePulses = FilterDatabase(PlottingDoublePulses, float(0.), IndexTable['CrossPolarization']) #removes fields where Polarization1 != Polarization2
  PlottingDoublePulses = FilterDatabase(PlottingDoublePulses, PulseDuration, IndexTable['pulses2']) # Removing fields were pulseDuration1 != pulseDuration2
  
  #### 2D PLOTTING Density(Color1, Color2) for one given set of laser parameters. We produce one plot like this for equal fields. 
  
  ## Selecting fields which have a strength close to a number indicated in input. 
  PlottingDoublePulses = FilterDatabaseGreaterThan(PlottingDoublePulses, FieldStrength*(1.-Tolerance), IndexTable['Efield1'])
  PlottingDoublePulses = FilterDatabaseLowerThan(PlottingDoublePulses, FieldStrength*(1.+Tolerance), IndexTable['Efield1'])
  PlottingDoublePulses = FilterDatabaseGreaterThan(PlottingDoublePulses, FieldStrength*(1.-Tolerance), IndexTable['Efield2'])
  PlottingDoublePulses = FilterDatabaseLowerThan(PlottingDoublePulses, FieldStrength*(1.+Tolerance), IndexTable['Efield2'])
  
  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq, wavelength2_uniq, CEP2_uniq, Efield2_uniq, pulses2_uniq= CountTheNumberOfLaserParameters(PlottingDoublePulses)
  
  print Header+"Total number of datasets: "+str(len(set(wavelength_uniq))*len(set(wavelength2_uniq)))
  DensityMax2D = np.empty((0,3)) #preparing the 2D plot
  for Color1 in wavelength_uniq: #selecting color1
    PlottingDoublePulses_SelectedColor1 = FilterDatabase(PlottingDoublePulses, Color1, IndexTable['wavelengths1'])
    for Color2 in wavelength2_uniq: #selecting color2
      # Selecting Color 2
      PlottingDoublePulses_SelectedColor2 = FilterDatabase(PlottingDoublePulses_SelectedColor1, Color2, IndexTable['wavelengths2'])
      
      # Constructing the density table
      runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingDoublePulses_SelectedColor2)
      
      print "Color1: "+str(Color1)+" m, Color2: "+str(Color2)+" m."
      if(len(Ne_final) >= 1):
        ToBeAdded=np.array([wavelengths1[0], wavelengths2[0], Ne_final[0]]) #BUG: removes some results, not clean enough!
        print ToBeAdded
        DensityMax2D=np.vstack((DensityMax2D, ToBeAdded))
        del ToBeAdded
        
  # Stack is built. Export, show and plot on a heat map. 
  print "** Info: exporting the sub-plot2D: Ne_final(wavelength, Field)..."
  print "** Info: useful to visualize available data with GNUplot. "
  filename = "Bicolor-Field"+str(round(FieldStrength*1E-9))+"VperNm-NeFinal"
  
  # Sorting the data before storing and plotting
  DensityMax2Dsorted = DensityMax2D[np.lexsort((DensityMax2D[:,2], DensityMax2D[:,1], DensityMax2D[:,0]))]
  #DensityMax2D.sort(axis=1) ##WARNING: THIS IS OUT OF CONTROL. DO NOT USE. 
  
  ExportToTxt(DensityMax2Dsorted, filename+".csv") #this is voluntarily unsorted to avoid damage
  
  #print "Sorted matrix?"
  #print DensityMax2Dsorted
  BrutX = 1E9*DensityMax2Dsorted[:,0]
  BrutY = 1E9*DensityMax2Dsorted[:,1]
  BrutZ = DensityMax2Dsorted[:,2]
  print len(BrutX)
  print len(BrutY)
  print len(BrutZ)
  
  X=list(set(BrutX))
  Y=list(BrutY) ##BUG: we have degeneracy, investigate why and remove the problem. 
  Z=DensityMax2Dsorted
  #print X
  #print Y
  
  #print len(X)
  #print len(Y)
  YY, XX = np.meshgrid(X, Y)
  #print np.shape(XX)
  #print np.shape(YY)
  #print np.shape(DensityMax2Dsorted)
 
  #print DensityMax2Dsorted
  #print np.shape(XX)
  #DensityMax2DsortedShape = np.reshape(DensityMax2Dsorted, np.shape(XX))
  
  # Interpolation to smooth the picture and colorize it more softly
  # BUG: function interp2d does not work! 
  #MeshPrecision=50
  #DensityMax2D_interp = interp2d(BrutX, BrutY, BrutZ)
  #Xi = np.arange(np.min(BrutX), np.max(BrutX), (np.max(BrutX)-np.min(BrutX))/MeshPrecision)
  #Yi = np.arange(np.min(BrutY), np.max(BrutY), (np.max(BrutY)-np.min(BrutY))/MeshPrecision)
  #YYi, XXi = np.meshgrid(Xi, Yi)
  ##print XXi
  #BicolorNeFinal=DensityMax2D_interp(Xi, Yi)
  
  # If the size is big enough, let's make a plot.
  if(len(DensityMax2D) >= 2):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.title(r'$N_{exc}(\lambda_1,\lambda_2)$, $E_{1,2}=$'+str(FieldStrength*1E-9)+' V/nm')
    plt.xlabel(r'Wavelength $\lambda_1$ (nm)')
    plt.ylabel(r'Wavelength $\lambda_2$ (nm)')
    #CS = plt.imshow(XX, YY, DensityMax2Dsorted) #, cmap=plt.cm.Blues)
    #CS = plt.contourf(XX, YY, DensityMax2Dsorted, cmap=plt.cm.Blues)
    #CS = ax.scatter(BrutX, BrutY, BrutZ) #, cmap=plt.cm.Blues) #WORKS
    #CS = ax.plot_wireframe(XX, YY, DensityMax2Dsorted) #, cmap=plt.cm.Blues) #FAILS: WHY do we have ground values ?! 
    #CS = ax.plot_surface(XX, YY, DensityMax2DsortedShape) #, cmap=plt.cm.Blues) #FAILS: WHY do we have ground values ?! 
    CS = ax.plot_trisurf(BrutX, BrutY, BrutZ, cmap=plt.cm.Blues) #WORKS
    
    #print np.shape(BrutX), np.shape(BrutY), np.shape(BrutZ)
    #print np.shape(Xi), np.shape(Yi)
    #print np.shape(XXi), np.shape(YYi)
    #print np.shape(BicolorNeFinal)
    
    #CS = ax.scatter(XXi, YYi, BicolorNeFinal) #, cmap=plt.cm.Blues) #WORKS
    
    #CS = ax.plot_trisurf(Xi, Yi, BicolorNeFinal, cmap=plt.cm.Spectral)
    
    #CS = ax.contour(XX, YY, Z, cmap=plt.cm.Blues) #WORKS
    #plt.colorbar(CS)
  
    plt.savefig(filename+'.eps')
    plt.savefig(filename+'.png')
    if(ShowFinalPlot==1):
      plt.show()
  
  del DensityMax2D, DensityMax2Dsorted
#}}}

# VALID but some degeneracies !
#plotEffectOfIntensity(PlottingTable, 10e-15, 0) #validated
plotEffectOfIntensity(PlottingTable, 20e-15, 0) #BUG: difference may originate from RisingTime !!
#plotEffectOfIntensity(PlottingTable, 30e-15, 0) #validated
#plotEffectOfIntensity(PlottingTable, 40e-15, 0) #validated
#plotEffectOfIntensity(PlottingTable, 50e-15, 0) #validated
plotEffectOfIntensity_PulseDuration(PlottingTable, 800e-9) #BUG: still a problem of shifted results
#plotEffectOfWavelength_Tolerance(PlottingTable) #considers a tolerance to plot various fields
plotEffectOfWavelength_2dMap(PlottingTable) #VALID works, but not enough data for now. 

# Developing...
#plotEffectOfCEP(PlottingTable)
# plotEffectOfFrozenXC(PlottingTable)
plotEffectOfDoublePulse_Color1_Color2_SingleFieldStrength(PlottingTable, 5E9) #VALID !!!! PRL !!!
#plotEffectOfDoublePulse_Efield1_Efield2_batchAllColors(PlottingTable) #BUGGY
#plotEffectOfDoublePulse_sameAmplitude_samePulseDurations(PlottingTable)
#plotEffectOfDoublePulse_diffWavelengths(PlottingTable)

print Header+"End of program."
exit()

