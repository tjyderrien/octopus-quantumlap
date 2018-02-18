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
# from scipy import interpolate
import matplotlib.pyplot as pplt
# from pylab import *
import sys
import decimal


# import libKeldysh as Keldysh
#from libAtomicUnits import *
from libPlottingSettings import *
from libPlotting import *
from libKeldysh import *
from libDatabase import *
from libExcitationDatabase import *

Header="[plotFinalQuantities] "

ShowFinalPlot              = 0

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

PlottingTable = GetExcitationDatabase(path1)
if(len(PlottingTable) > 1): 
  print Header+"** Info: capturing data from old PATH1."
else:
  PlottingTable = GetExcitationDatabase(path2)
  print Header+"** Info: capturing data from novel path PATH2."

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
SizeX = 8.
SizeY = 17.

pplt.figure(figsize=(SizeX,SizeY))

## We change strategy: now we plot all the available data based on the number of parameters we have detected. The input parameters will not define what we have to plot anymore. We prepare all the plots, which will simplify the understanding of this file. 
def plotEffectOfIntensity(PlottingTable):
  print "#====== Plot: SINGLE PULSE, SAME WAVELENGTHS, SAME POLARIZATIONS: final density as function of field ======"
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

  ## prints the current state of selection
  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq = CountTheNumberOfLaserParameters(PlottingSinglePulses)

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
    PlottingSinglePulses_wavelength = FilterDatabase(PlottingSinglePulses, 800e-9, IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts.

    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingSinglePulses_wavelength)

    pplt.plot(FieldToIntensity(Efield1), Ne_final, 'or', label=r'$N_{exc}$, $\lambda=$'+str(round(wavelength_value*1E9))+' nm')
  #}}}
  ax1.legend(loc=2, borderaxespad=0.) #bbox_to_anchor=(1.05, 1), 
  #pplt.xlabel(r'Peak intensity (W/m$^2$)')
  pplt.ylabel(r'Final exc. density (m$^{-3}$)')
#  pplt.ylim((1E16,4.*5e28))
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  
  ax12 = ax1.twinx()
  for wavelength_value in wavelength_uniq: #{{{
    # select the corresponding wavelength in results
    PlottingSinglePulses_wavelength = FilterDatabase(PlottingSinglePulses, 800e-9, IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts.

    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, Energy_max, W_pi, PulseType = ExtractData(PlottingSinglePulses_wavelength)

    pplt.plot(FieldToIntensity(Efield1), Ne_final/pulses1, 'ro') #, label=r'$N_{exc}/\tau$, $\lambda=$'+str(round(wavelength_value*1E9))+' nm')
  #}}}
  ax12.set_yscale('log')
  pplt.ylabel('$w_{PI}$ (m$^{-3}$.s$^{-1}$)')
  
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
  ax12.legend(loc=4, borderaxespad=0.)
  ax22.legend(loc=1, borderaxespad=0.)
  
  # Saving the plot
  ExportName = "EffectOfWavelength"
  pplt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.eps')
  pplt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.png')
  print Header+"** Saved: Metaplot-SinglePulse-SingleColor-EffectOfWavelength.png."
  #if(ShowFinalPlot==1):
    #pplt.show()

  del ax1, ax12, ax2, ax22, PlottingSinglePulses, wavelength_uniq, CEP1_uniq
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

  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq = CountTheNumberOfLaserParameters(PlottingSinglePulses)
  
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

  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq = CountTheNumberOfLaserParameters(PlottingSinglePulses)

  #if(len(CEP1_uniq) == 1): #{{{ #single color simulations
  # Complexity level: 2


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

##TODO: shall we make a multi-wavelength figure?
##TODO: - If wavelength2!=wavelength1, yes.
##TODO: - If len(set(wavelength1)) > 1, yes.  

#del PlottingTable

def EffectOfWavelength(PlottingTable):
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
  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq = CountTheNumberOfLaserParameters(PlottingSinglePulses)
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
  pplt.savefig('Metaplot-SinglePulse-SingleColor-EffectOfWavelength.eps')
  pplt.savefig('Metaplot-SinglePulse-SingleColor-EffectOfWavelength.png')
  #if(ShowFinalPlot==1):
    #pplt.show()

  del ax1, PlottingSinglePulses
  #}}}

def EffectOfDoublePulse():
  print "#====== Plot: DOUBLE PULSE, SAME WAVELENGTH, SAME POLARIZATION: final density as function of delay and field ======"
  print "Jammed. Have to reorganize the double-pulse plotting."
  #if(EnableDoublePulseAnalysis==1): #{{{
    #if(wavelength_numbers == 1): 
      ##Plot 3: Final density as function of delay
      #pplt.title('Density as function of delay')
      #pplt.scatter(Delay*1E15, Ne_final, s=5.+(FieldToIntensity(Efield1)/1E14)*10., c=10*np.rint(CrossPolarization*180./pi), label=r'$N_{exc}$[Delay], Pola='+str(CrossPolarization[0]*180./pi)+"deg")
      ## ax3.legend(bbox_to_anchor=(0.00, -0.1), loc=3, borderaxespad=0.)
      #pplt.xlabel('Double pulse delay (fs)')
      #pplt.ylabel(r'Final density (m$^{-3}$)')
      ## plt.ylim((1E16,4.*5e28))
      ## plt.xlim((0.5E-1, 35.))
      ## ax3.set_xscale("log")
      ## ax3.set_yscale("log")

      ##Plot 4: Irreversibility as function of delay
      #ax4 = plt.subplot(plotD)
      #pplt.title('Irreversibility factor')
      #pplt.scatter(Delay*1E15, irreversibility, s=5.+(FieldToIntensity(Efield1)/1E14)*10., c=10*np.rint(CrossPolarization*180./pi), label='Irreversibility[Delay]')
      ## ax4.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
      #pplt.xlabel('Double pulse delay (fs)')
      #pplt.ylabel('Irreversibility factor')
      #pplt.ylim((0,1))
    ##}}}
  ## Final saving
  #figmode=''

  #if(EnableDoublePulseAnalysis == 1):
    #figmode = figmode+'-DoublePulse'
  #if(EnableInducedFieldAnalysis == 1):
    #figmode = figmode+'-InducedFields'

  #figname = 'Metaplot'+figmode
  #plt.savefig(figname+'.eps', bbox_inches='tight')
  #plt.savefig(figname+'.png', bbox_inches='tight')
  #if(ShowFinalPlot==1):
    #plt.show()

  ## Plot as function of wavelength for one given peak intensity
  ## Plot as function of peak intensity for one given wavelength.
  ## Plot as function of double pulse delay

  #print "====== Plot: DOUBLE PULSE, SAME WAVELENGTH, SAME POLARIZATION: final density as function of delay and field ======"
  ##Filter: double pulses, same wavelengths, same polarization
  #PlottingDoublePulses = PlottingTable #get the full database
  #PlottingDoublePulses = FilterDatabaseGreaterThan(PlottingDoublePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak < 1 V/m. #works
  #PlottingDoublePulses = FilterDatabase(PlottingDoublePulses, wavelengths1, IndexTable['wavelengths2']) #removes fields where wavelength2 != wavelength1
  #PlottingDoublePulses = FilterDatabase(PlottingDoublePulses, float(0.), IndexTable['CrossPolarization']) #removes fields where wavelength2 != wavelength1
  #runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Ne_max, Ne_final, irreversibility, Energy_max = ExtractData(PlottingDoublePulses)

  ##ax1 = plt.subplot(plotA)

  ## Number of sweeped parameters: 2
  #plt.title('Final density (peak intensity, delay)')
  #print "** 3d plot: Delay*       size = "+str(len(Delay_uniq))
  #print "** 3d plot: wavelength1* size = "+str(len(wavelengths1_uniq))
  ##  Ne_final_3d=np.reshape(Ne_final, (len(Efield1_uniq), len(Delay_uniq)))
  ## Ne_final_3d_interp = interpolate.interp2d(Efield1, Delay, Ne_final)
  ## Efield1_interp = np.arange(1E8, 1E10, 1E8)
  ## Delay_interp   = np.arange(0e0, 50e-15, 10e-15)
  ## Ne_final_3d = Ne_final_3d_interp(Efield1_interp, Delay_interp) 
  ## plt.contourf(Efield1_interp, Delay_interp, Ne_final_3d, tri=True)
  #pplt.ylabel('Delay (fs)')
  #pplt.xlabel('Field ampl. (V/m)')
  ##pplt.scatter(FieldToIntensity(Efield1), Delay*1E15, s=34+30*Ne_final/Density, c=Ne_final/Density)
  #pplt.scatter(FieldToIntensity(Efield1), Delay*1E15, c=Ne_final/Density)
  #pplt.savefig('Metaplot-DoublePulse.eps')
  #pplt.savefig('Metaplot-DoublePulse.png')
  #pplt.show()

plotEffectOfIntensity(PlottingTable) #validated
#plotEffectOfCEP(PlottingTable)
# plotEffectOfFrozenXC(PlottingTable)

print Header+"End of program."
exit()

