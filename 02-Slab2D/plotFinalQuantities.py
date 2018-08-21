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
import itertools

import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.legend_handler import HandlerLine2D
import matplotlib.cm as cm
import mpl_toolkits.mplot3d as mplt 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as tri
#from matplotlib.tri import Triangulation, UniformTriRefiner
from scipy.interpolate import interp2d
from scipy.optimize import curve_fit
from colorama import Fore
from colorama import Style

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

ShowFinalPlot              = 0

Density=5E28 #specific to Si #TODO: compute automatically? 

#Egap_SI = 2.56*e #LDA: 2.56 eV, TB09: 3.04 eV. [Waroquiers and Gonze, Phys. Rev. B 87, 075121 (2013)] #TODO: this will depend on which functional we have used! TB09 or LDA ! 

# Capture the data
# For histocal reasons, two locations are possilbe.
path1='*/*/output_iter/' #valid for the runs computed from a same k-grid folder. 
path2='*/*/td.general/'

PlottingTable = libExcitationDatabase.GetExcitationDatabase(path2) #attempts to write a new DB. 
#if(len(PlottingTable) > 1): 
  #print Header+"** Info: capturing data from old PATH1." #BUG: destroys the previous database! 
#else:
  #PlottingTable = GetExcitationDatabase(path2)
  #print Header+"** Info: capturing data from novel path PATH2."

print PlottingTable

PlottingTable = FilterDatabase(PlottingTable, 0, IndexTable['FieldInduced']) #removes all induced fields data

runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingTable)

print Header+"** Preparing the summary plots..."

#### =============== PLOTTING STARTS HERE ==================
# The filtering will be done on the fly. Then we are sure of 
# what we plot in each figure. 

### Defining plotting structure

# Observables:        {N_exc^max, N_exc^final, irreversibility, TotalCurrent_final, Energy_max, Energy_final}

#### Plotting: parameters for controlling the canvas
# if(EnableDoublePulseAnalysis==0):
  # 2 plots
plotA=211; plotB=212; #plotC=313
# elif (wavelength_numbers > 1):
# plotA=311; plotB=312; plotC=313;
# else: #if(EnableDoublePulseAnalysis==1):
# 3 plots
# plotA=311; plotB=312; plotC=313;
# 4 plots 
# plotA=221; plotB=222; plotC=223; plotD=224;

# Size of main canvas
GoldNumber = ((1.+np.sqrt(5.))/2.)
SizeX = 8.
SizeY = SizeX / GoldNumber

#pplt.figure(figsize=(SizeX,SizeY))

## Plot Nexc(t, single tau, 1 line per wavelength)
def plotEffectOfRiseTime(PlottingTable, pulseWavelength=800e-9, pulseDuration=20e-15, Efield1=5E9, PulseShape=3, ShowFinalPlot=0, Egap_SI=2.56*e):
  print "{Fore.GREEN}#====== Plot: SINGLE PULSE, EFFECT OF THE RISE TIME, 1 WAVELENGTH, 1 POLARIZATIONS: final density as function of RISETIME ======{Style.RESET_ALL}"
  PlottingSinglePulses = PlottingTable #get the full database
  PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
  PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, pulseWavelength, IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts. #TODO: this field should be making a loop in available wavelengths given in wavelengths1_uniq.
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0, IndexTable['FreezeHXC']) #keeps only non-frozen XC. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, PulseShape, IndexTable['PulseType']) #keeps only non-frozen XC. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, pulseDuration, IndexTable['pulses1']) #removes fields where tau != pulseDuration

  PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, Efield1*1.01, IndexTable['Efield1'])
  PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, Efield1*0.99, IndexTable['Efield1'])

  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = CountTheNumberOfLaserParameters(PlottingSinglePulses)
  
  pplt.figure(figsize=(SizeX,SizeY))
  ax1 = pplt.subplot(111) #111 or plotA if using plotB
  pplt.title(r'$n_{exc}(RiseCycles)$, $\tau=$'+str(int(pulseDuration*1E15))+' fs.')
  # Number of sweeped parameters: 3 (Efield, Delay, CrossPolarization)
  #pplt.scatter(FieldToIntensity(Efield1), Ne_final, c=Delay, s=(30.+CrossPolarization*50.), label=r'$N_{exc}$[Field]')
  plot1 = [] #np.empty((0,0))
  for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
    PlottingSinglePulses_CEP = FilterDatabase(PlottingSinglePulses, CEP1_value, IndexTable['CEP1']) #select the local CEP
    #TODO: filter by field strength ! Maybe take a piece of code from tolerance-based selector in another function. #TODO: Use subfunctions otherwise we will get crazy. 
    ## prints the current state of selection and stores the values of sweeped parameters
    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_CEP)
    
    #plot_field, = ax1.plot(FieldToIntensity(Efield1), Ne_final, '-', label=r'$N_{exc}^{final}$, $\phi_{CE}$='+str(round(CEP1_value*180./pi)))
    #lists = sorted(itertools.izip(*[FieldToIntensity(Efield1), Ne_final]))
    lists = sorted(itertools.izip(*[RiseTime, Ne_max, Ne_final])) #NOTE: Ne_max makes more sense than Ne_final
    RiseTime_sorted, Ne_max_sorted, Ne_final_sorted = list(itertools.izip(*lists))
    #print Header+"Ready to plot: "
    #print "** Info: Sorting data..."
    #print Efield1_sorted, Ne_final_sorted
    plot_field, = ax1.plot(RiseTime_sorted, Ne_final_sorted, '-', label=r'$N_{exc}^{final}$, $\phi_{CE}$='+str(round(CEP1_value*180./pi)))
    print Header+"plotEffectOfIntensity_CEP(): adding CEP "+str(round(CEP1_value*180./pi))+" to the the plot."
    print plot_field
    #print Header+"** Debug info: "+str(np.shape(plot_field))
    plot1 += [plot_field] # does cumulate each CEP now. Type: Line2D, not np.array. 
    
    #ax1.legend(loc=4, borderaxespad=0.) #bbox_to_anchor=(1.05, 1), 
  labels1 = [l.get_label() for l in plot1]
  ax1.legend(plot1, labels1, loc='best', fancybox=True, framealpha=0)
  plt.xlabel(r'Rise cycles')
  plt.ylabel(r'Maximum density (m$^{-3}$)')
  #  pplt.ylim((1E16,4.*5e28))
  #ax1.set_xscale('log')
  #ax1.set_yscale('log')
  #plot1_list = plot1.tolist()
    #print Efield1, CEP1_value
    
  #Adding axis of gammaKeldysh, but for all CEPs !!
  #runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses)
  
  ax12 = ax1.twinx()
  for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
    PlottingSinglePulses_CEP = FilterDatabase(PlottingSinglePulses, CEP1_value, IndexTable['CEP1']) #select the local CEP
    
    ## prints the current state of selection and stores the values of sweeped parameters
    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_CEP)
    
    #plot_field, = ax1.plot(FieldToIntensity(Efield1), Ne_final, '-', label=r'$N_{exc}^{final}$, $\phi_{CE}$='+str(round(CEP1_value*180./pi)))
    #lists = sorted(itertools.izip(*[FieldToIntensity(Efield1), Ne_final]))
    lists = sorted(itertools.izip(*[RiseTime, Ne_max, Ne_final, pulses1])) #NOTE: Ne_max makes more sense than Ne_final
    RiseTime_sorted, Ne_max_sorted, Ne_final_sorted, pulses1_sorted = list(itertools.izip(*lists))
    plt.plot(RiseTime_sorted, np.divide(Ne_final_sorted, pulses1_sorted), 'o') #, label=r'$N_{exc}/\tau$') #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
  
  #ax12.set_yscale('log')
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
  #ax1.legend(bbox_to_anchor=(1.15, 1), loc=2, borderaxespad=0.)
  #ax22.legend(loc=1, borderaxespad=0.)

  # Saving the plot
  ExportName = "EffectOfRiseTime"
  plt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.eps')
  plt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.png')
  if(ShowFinalPlot==1):
    pplt.show()

  del ax1, PlottingSinglePulses, wavelength_uniq, wavelengths1
#}}}

## Plot Nexc(t, single tau, 1 line per wavelength)
def plotEffectOfIntensity_wavelengths(PlottingTable, pulseDuration=20e-15, CEP=0., RestrictWavelength1=False, wavelength_value=800e-9, PulseType=3, RiseTime=2, ShowFinalPlot=0, Egap_SI = 2.56*e):
  PlotVirtualExcitation = False
  print "{Fore.GREEN}#====== Plot: SINGLE PULSE, VARIOUS WAVELENGTHS, SAME PULSE DURATION, SAME POLARIZATIONS: final density as function of field ====== {Style.RESET_ALL}"
  #Filter: double pulses, same wavelengths, same polarization
  PlottingSinglePulses = PlottingTable #get the full database
  PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
  PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, CEP, IndexTable['CEP1']) #removes fields where CEP != 0. #NOTE: Can be done for Si-LDA. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['FreezeHXC']) #removes fields where FreezeHXC != 0. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['FieldInduced']) #removes fields where FreezeHXC != 0. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 1., IndexTable['TheoryLevel']) #removes runs where TheoryLevel != dft. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, PulseType, IndexTable['PulseType']) #keeps only pulses of type StH (3)
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, RiseTime, IndexTable['RiseTime']) #keeps only pulses of type STH RiseTime
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, pulseDuration, IndexTable['pulses1']) #keeps only pulses of selected duration in argument
  if(RiseTime >= 1): #select by RiseTime only if it is defined
      PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, RiseTime, IndexTable['RiseTime'])
  
  ## prints the current state of selection
  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq, wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = CountTheNumberOfLaserParameters(PlottingSinglePulses)
  if(PlotVirtualExcitation==True): 
      plotA=211; plotB=212
  else: 
      plotA=111; plotB=111
  # Simple pulse, single color
  pplt.figure()
  #if(len(wavelength_uniq) == 1): #{{{ #single color simulations
  # Complexity level: 1
  
  ax1 = pplt.subplot(plotA)
  
  pplt.title('Density as function of peak intensity')
  # Number of sweeped parameters: 3 (Efield, Delay, CrossPolarization)
  #pplt.scatter(FieldToIntensity(Efield1), Ne_final, c=Delay, s=(30.+CrossPolarization*50.), label=r'$N_{exc}$[Field]')
  for wavelength_value in wavelength_uniq: #{{{ #TODO: wavelength_uniq should be sorted (visible in the caption)
    #if(RestrictWavelength1==False): 
      # select the corresponding wavelength in results
    PlottingSinglePulses_wavelength = FilterDatabase(PlottingSinglePulses, wavelength_value, IndexTable['wavelengths1']) #keeps all wavelengths wavelenghts.
    #else: #selection only 1 wavelength #BUG: DOES NOT WORK, and even attributes to other wavelength what was done for 800 nm
      #PlottingSinglePulses_wavelength = FilterDatabase(PlottingSinglePulses, 800e-9, IndexTable['wavelengths1']) #keeps only 800 nm

    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_wavelength)
    lists = sorted(itertools.izip(*[wavelengths1,FieldToIntensity(Efield1), Ne_final, pulses1]))
    wavelengths_sorted, Intensity_sorted, Ne_final_sorted, pulses1_sorted = list(itertools.izip(*lists))
    pplt.plot(Intensity_sorted, Ne_final_sorted, 'o', label=r'$N_{exc}$, $\lambda=$'+str(round(wavelength_value*1E9))+' nm')
  #}}}
  ax1.legend(loc='best', fancybox=True, framealpha=0) #bbox_to_anchor=(1.05, 1), 
  if(not PlotVirtualExcitation): 
    pplt.xlabel(r'Peak intensity (W/m$^2$)')
  pplt.ylabel(r'Final exc. density (m$^{-3}$)')
#  pplt.ylim((1E16,4.*5e28))
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  
  ax12 = ax1.twinx()
  for wavelength_value in wavelength_uniq: #{{{
    # select the corresponding wavelength in results
    if(RestrictWavelength1==False): 
      PlottingSinglePulses_wavelength = FilterDatabase(PlottingSinglePulses, 800e-9, IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts.
    else:
      PlottingSinglePulses_wavelength = FilterDatabase(PlottingSinglePulses, wavelength_value, IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts.

    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_wavelength)

    lists = sorted(itertools.izip(*[wavelengths1,FieldToIntensity(Efield1), Ne_final, pulses1]))
    wavelengths_sorted, Intensity_sorted, Ne_final_sorted, pulses1_sorted = list(itertools.izip(*lists))

    pplt.plot(Intensity_sorted, np.divide(Ne_final_sorted, (np.multiply(2.,pulses1_sorted))), 'o') #, label=r'$N_{exc}/\tau$, $\lambda=$'+str(round(wavelength_value*1E9))+' nm')
  #}}}
  ax12.set_yscale('log')
  pplt.ylabel('$w_{PI}$ (m$^{-3}$.s$^{-1}$)')
  
  # Second panel of the figure
  if(PlotVirtualExcitation):
    ax2 = pplt.subplot(plotB)
    pplt.title('rVirtual exc. $\zeta$')
    pplt.plot(FieldToIntensity(Efield1), irreversibility, 'ob', label='Virtual exc.')
    pplt.xlabel(r'Intensity (W/m$^2$)')
    pplt.ylabel(r'Virtual exc. $\zeta$')
    pplt.ylim((0,1))
    ax2.set_xscale('log')
    ax22 = ax2.twinx()
    for wavelength_value in wavelength_uniq:
      # select the corresponding wavelength in results
      PlottingSinglePulses_wavelength = FilterDatabase(PlottingSinglePulses, wavelength_value, IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts.

      runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_wavelength)
      pplt.plot(FieldToIntensity(Efield1), gammaKeldysh(Egap_SI*e, meff, Efield1, wavelengths1), 'or', label=r'Adiabaticity, $\lambda=$'+str(round(wavelength_value*1E9))+' nm')
    pplt.ylabel(r'Adiabadicity $\gamma$')
    #pplt.ylim((0,1))
    ax22.set_yscale('log')
      
    # ax2.set_yscale("log")
    ax2.legend(loc='best', fancybox=True, framealpha=0)
    # ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax22.legend(loc='best', fancybox=True, framealpha=0)
  ax12.legend(loc='best', fancybox=True, framealpha=0)
  
  # Saving the plot
  ExportName = "EffectOfIntensity-"+str(int(pulseDuration*1E15))+"fs"
  pplt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.eps')
  pplt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.png')
  print Header+"** Saved: Metaplot-SinglePulse-SingleColor-"+ExportName+".eps|png."
  if(ShowFinalPlot==1):
    pplt.show()

  del ax1, ax12, PlottingSinglePulses, wavelength_uniq, CEP1_uniq
  if(PlotVirtualExcitation): 
    del ax2, ax22
#}}}

## Plots N_exc(t; 1 plot per tau, single lambda)
def plotEffectOfIntensity_PulseDuration(PlottingTable, pulseWavelength=800e-9, CEP0=0., PulseType=3, RiseCycles=2, PlotWpiAxis=False, ShowFinalPlot=0, Egap_SI=2.56*e):
  PlotVirtualExcitation = True
  PlotEnergy = True
  EnableKeldysh = True
  meff = 0.2226
  print "#{Fore.GREEN}====== Plot: SINGLE PULSE, SAME WAVELENGTHS, SAME POLARIZATIONS: final density as function of field and pulse duration ======{Style.RESET_ALL}"
  #Filter: double pulses, same wavelengths, same polarization
  PlottingSinglePulses = PlottingTable #get the full database
  PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
  PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, CEP0, IndexTable['CEP1']) #removes fields where CEP != 0. #NOTE: Can be done for Si-LDA. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['FreezeHXC']) #removes fields where FreezeHXC != 0. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['FieldInduced']) #removes fields where FreezeHXC != 0. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 1., IndexTable['TheoryLevel']) #removes runs where TheoryLevel != dft. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, PulseType, IndexTable['PulseType']) #keeps only pulses of type StH (3)
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, CEP0, IndexTable['CEP1'])
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, RiseCycles, IndexTable['RiseTime']) #keeps only pulses of a selected RiseCycles
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, pulseWavelength, IndexTable['wavelengths1']) #keeps only pulses of type STH (2)

  ## prints the current state of selection
  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = CountTheNumberOfLaserParameters(PlottingSinglePulses)

  # Simple pulse, single color
  pplt.figure(figsize=(GoldNumber*SizeX,GoldNumber*SizeY))
  #if(len(wavelength_uniq) == 1): #{{{ #single color simulations
  if(PlotEnergy and PlotVirtualExcitation):
    #plotA = 311; plotB = 312; plotC = 313 #Horizontal
    plotA = 131; plotB = 132; plotC = 133 #Vertical
  elif(PlotVirtualExcitation):
    plotA = 211; plotB = 212
  else: 
    plotA = 111; plotB = 111
  
  EnergyMin = 1e-5 #-1E-7
  
  
  ax1 = pplt.subplot(plotA)
  pplt.title(r'$\lambda = $'+str(round(pulseWavelength*1E9))+' nm')
  
  plot1 = []
  for pulse_value in pulses_uniq: #{{{
    # select the corresponding wavelength in results
    PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses, pulse_value, IndexTable['pulses1'])
    PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses_pulses, pulse_value, IndexTable['pulses2']) #keeps only pulses of the same duration.

    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_pulses)
    
    print Header+"** Info: saving the data in CSV file for fitting. "
    filename = "EffectOfPulseIntensity-Duration-"+str(pulseWavelength*1E9)+"nm-tau"+str(int(pulse_value*1E15))+"fs"
    print Header+"** Info: saving in "+filename+".csv"
    ExportToTxt(PlottingSinglePulses_pulses, filename+".csv")
    lists = sorted(itertools.izip(*[pulses1, Efield1, FieldToIntensity(Efield1), Ne_final]))
    pulses1_sorted, Efield1_sorted, Intensity_sorted, Ne_final_sorted = list(itertools.izip(*lists))
    #print Header+"Ready to plot: "
    
    plot_field, = ax1.plot(Intensity_sorted, Ne_final_sorted, '^', label=r'$n_{exc}^{tddft}$, $\tau_1=$'+str(int(pulse_value*1E15))+' fs')
    plot1 += [plot_field]
    
    #==== CURVE FITTING SHOULD COME HERE === as it is for each pulse duration. 
    print Header+"** Info: Computation of multi-photon excitation"
    
    ### Multiphotonic absorption law (used for fitting td-dft results)
    #def MultiPhotonAbsorption(intensity, sigma_n):
        #n_photons=2 # is defined by inheriting its value
        #omega = 2.*np.pi*c/wavelengths1[0] #how come this value is taken outside, but not my n_photon quantity?
        ##n_photons #gets the number of photons in upper levels
        #return sigma_n * intensity**n_photons / (n_photons*hbar*omega) * pulses1_sorted[0]

    ## 1-photon absorption fitting function
    def OnePhotonAbsorption(intensity, sigma_n): 
    #{{{
      n=1; #wavelengths1=800e-9; #is captured by parental legacy
      omega = 2.*np.pi*c/wavelengths1[0]
      #nexc_mpi = MultiPhotonAbsorption(intensity, sigma_n, n, wavelengths1[0], pulses1_sorted[0])
      nexc_mpi = sigma_n * intensity**n / (n*hbar*omega) * pulses1_sorted[0]
      return nexc_mpi
    #}}}
    
    ## 2-photon absorption fitting function
    def TwoPhotonAbsorption(intensity, sigma_n): 
    #{{{
      n=2; #wavelengths1=800e-9; #is captured by parental legacy
      omega = 2.*np.pi*c/wavelengths1[0]
      #nexc_mpi = MultiPhotonAbsorption(intensity, sigma_n, n, wavelengths1[0], pulses1_sorted[0])
      nexc_mpi = sigma_n * intensity**n / (n*hbar*omega) * pulses1_sorted[0]
      return nexc_mpi
    #}}}
    
    ## 3-photon absorption fitting function
    def ThreePhotonAbsorption(intensity, sigma_n): 
    #{{{
      n=3; #wavelength=800e-9; tau_FWHM=20e-15; 
      omega = 2.*np.pi*c/wavelengths1[0]
      tau_FWHM = pulses1_sorted[0]
      #nexc_mpi = MultiPhotonAbsorption(intensity, sigma_n, n, wavelength, tau_FWHM)
      nexc_mpi = sigma_n * intensity**n / (n*hbar*omega) * tau_FWHM
      return nexc_mpi
    #}}}
    
    ## 2-photon and 3-photon absorption fitting function
    def TwoAndThreePhotonAbsorption(intensity, sigma_2, sigma_3): 
    #{{{
      n=2; #wavelengths1=800e-9; #is captured by parental legacy
      omega = 2.*np.pi*c/wavelengths1[0]
      #nexc_mpi = MultiPhotonAbsorption(intensity, sigma_n, n, wavelengths1[0], pulses1_sorted[0])
      nexc_mpi = sigma_2 * intensity**2 / (2*hbar*omega) * pulses1_sorted[0] + sigma_3 * intensity**3 / (3*hbar*omega) * pulses1_sorted[0]
      return nexc_mpi
    #}}}
    
    ## Filters a part of the data and fit using the given function
    def FitTwoPhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, n_photons, TwoPhotonIntensityRegime_min, TwoPhotonIntensityRegime_max): #TODO: pass function name in argument
    #{{{
        print Header+"** Info: Filtering the fitting region..."
        Intensity_sorted_np = np.array(Intensity_sorted) #filtering works for NumPy array, not for lists. 
        #condition = (Intensity_sorted_np > TwoPhotonIntensityRegime_min) and (Intensity_sorted_np < TwoPhotonIntensityRegime_max)
        Intensity_sorted_mask = Intensity_sorted_np[(Intensity_sorted_np > TwoPhotonIntensityRegime_min) & (Intensity_sorted_np < TwoPhotonIntensityRegime_max)]
        
        Ne_final_sorted_np  = np.array(Ne_final_sorted)
        Ne_final_sorted_mask  = Ne_final_sorted_np[(Intensity_sorted_np > TwoPhotonIntensityRegime_min) & (Intensity_sorted_np < TwoPhotonIntensityRegime_max)]
        print Intensity_sorted_mask, Ne_final_sorted_mask #Works well
        
        TwoPhotonAbsorption_params = curve_fit(TwoPhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
        #TwoPhotonAbsorption_params = curve_fit(MultiPhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
        print Header+"** Ionization parameters for tau="+str(int(1e15*pulse_value))+" fs: "+str(TwoPhotonAbsorption_params[0])
        
        print Header+"** Info: adding the fitted curve to the plot"
        Nexc_2PI = TwoPhotonAbsorption(Intensity_sorted_mask, TwoPhotonAbsorption_params[0])
        #print Nexc_2PI
        plot_field, = ax1.plot(Intensity_sorted_mask, Nexc_2PI, '--', label=r'Fit $\sigma_2=$'+str('%.4g'%TwoPhotonAbsorption_params[0][0])+r', $\tau_1=$'+str(int(pulse_value*1E15))+' fs') #NOTE add value of sigma_2 here
        del Nexc_2PI, TwoPhotonAbsorption_params, Intensity_sorted_mask
        return plot_field
    #}}}
    
        ## Filters a part of the data and fit using the given function
    def FitOnePhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, n_photons, OnePhotonIntensityRegime_min, OnePhotonIntensityRegime_max): #TODO: pass function name in argument
    #{{{
        print Header+"** Info: Filtering the fitting region..."
        Intensity_sorted_np = np.array(Intensity_sorted) #filtering works for NumPy array, not for lists. 
        #condition = (Intensity_sorted_np > OnePhotonIntensityRegime_min) and (Intensity_sorted_np < OnePhotonIntensityRegime_max)
        Intensity_sorted_mask = Intensity_sorted_np[(Intensity_sorted_np > OnePhotonIntensityRegime_min) & (Intensity_sorted_np < OnePhotonIntensityRegime_max)]
        
        Ne_final_sorted_np  = np.array(Ne_final_sorted)
        Ne_final_sorted_mask  = Ne_final_sorted_np[(Intensity_sorted_np > OnePhotonIntensityRegime_min) & (Intensity_sorted_np < OnePhotonIntensityRegime_max)]
        print Intensity_sorted_mask, Ne_final_sorted_mask #Works well
        
        OnePhotonAbsorption_params = curve_fit(OnePhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
        #TwoPhotonAbsorption_params = curve_fit(MultiPhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
        print Header+"** Ionization parameters for tau="+str(int(1e15*pulse_value))+" fs: "+str(OnePhotonAbsorption_params[0])
        
        print Header+"** Info: adding the fitted curve to the plot"
        Nexc_1PI = OnePhotonAbsorption(Intensity_sorted_mask, OnePhotonAbsorption_params[0])
        #print Nexc_2PI
        plot_field, = ax1.plot(Intensity_sorted_mask, Nexc_1PI, '--', label=r'Fit $\sigma_1=$'+str('%.4g'%OnePhotonAbsorption_params[0][0])+r', $\tau_1=$'+str(int(pulse_value*1E15))+'fs') #NOTE add value of sigma_1 here
        del Nexc_1PI, OnePhotonAbsorption_params, Intensity_sorted_mask
        return plot_field
    #}}}
    
    ## ==== FITTING THE 3-photon ionization plot
    def FitThreePhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, n_photons, ThreePhotonIntensityRegime_min, ThreePhotonIntensityRegime_max): #TODO: pass function name in argument instead of copying this function... :-(
    #{{{
        print Header+"** Info: Filtering the fitting region..."
        Intensity_sorted_np = np.array(Intensity_sorted) #filtering works for NumPy array, not for lists. 
        #condition = (Intensity_sorted_np > TwoPhotonIntensityRegime_min) and (Intensity_sorted_np < TwoPhotonIntensityRegime_max)
        Intensity_sorted_mask = Intensity_sorted_np[(Intensity_sorted_np > ThreePhotonIntensityRegime_min) & (Intensity_sorted_np < ThreePhotonIntensityRegime_max)]
        
        Ne_final_sorted_np  = np.array(Ne_final_sorted)
        Ne_final_sorted_mask  = Ne_final_sorted_np[(Intensity_sorted_np > ThreePhotonIntensityRegime_min) & (Intensity_sorted_np < ThreePhotonIntensityRegime_max)]
        print Intensity_sorted_mask, Ne_final_sorted_mask #Works well
        
        ThreePhotonAbsorption_params = curve_fit(ThreePhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
        #TwoPhotonAbsorption_params = curve_fit(MultiPhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
        print Header+"** Ionization parameters for tau="+str(int(1e15*pulse_value))+" fs: "+str(ThreePhotonAbsorption_params[0])
        
        print Header+"** Info: adding the fitted curve to the plot"
        Nexc_3PI = ThreePhotonAbsorption(Intensity_sorted_mask, ThreePhotonAbsorption_params[0])
        #print Nexc_2PI
        plot_field, = ax1.plot(Intensity_sorted_mask, Nexc_3PI, '--', label=r'Fit $\sigma_3=$'+str(ThreePhotonAbsorption_params[0][0])+r', $\tau_1=$'+str(int(pulse_value*1E15))+'fs') #NOTE adds value of sigma_3 here
        del Nexc_3PI, ThreePhotonAbsorption_params, Intensity_sorted_mask
        return plot_field
    #}}}
    
    ThreePhotonIntensityRegime_min = 0e0
    ThreePhotonIntensityRegime_max = 0e0 #W/m2
    
    TwoPhotonIntensityRegime_min = 0e0
    TwoPhotonIntensityRegime_max = 1.3e15 #W/m2
    
    OnePhotonIntensityRegime_min = 1e12
    OnePhotonIntensityRegime_max = 1.3e15 #W/m2
    
    omega = 2.*np.pi*c/wavelengths1[0]
    NumberOfPhotons = Egap_SI / (hbar * omega)
    NumberOfPhotons_round = round(NumberOfPhotons,0)
    print Header+"** INFO: NUMBER OF PHOTONS = "+str(NumberOfPhotons)+", rounded to: "+str(NumberOfPhotons_round)
    
    try:
        if(NumberOfPhotons_round >= 2 and NumberOfPhotons_round < 3): 
            print Header+"** Info: fitting with 2-photon absortion..."
            plot_field = FitTwoPhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, 2, TwoPhotonIntensityRegime_min, TwoPhotonIntensityRegime_max)
            plot1 += [plot_field]
        elif (NumberOfPhotons_round >= 1 and NumberOfPhotons_round < 2): 
            print Header+"** Info: fitting with 1-photon absortion..."
            plot_field = FitOnePhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, 1, OnePhotonIntensityRegime_min, OnePhotonIntensityRegime_max)
            plot1 += [plot_field]
        else: 
            print Header+"** Info: fitting with 3-photon absortion..."    
            plot_field = FitThreePhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, 3, ThreePhotonIntensityRegime_min, ThreePhotonIntensityRegime_max)
            plot1 += [plot_field]
    except: 
        print Header+"** Fit using N-photon absorption failed. "
    # NOTE: this is NOT 3-photon absorption regime, this is NON-PERTURBATIVE regime. Transition to tunneling. 
    
    print Header+"** Info: fitting the non-perturbative regime [DISABLED FOR NOW]..."
    #try: 
        #plot_field = FitTwoPhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, 2, ThreePhotonIntensityRegime_min, ThreePhotonIntensityRegime_max)
        #plot1 += [plot_field]
    #except: 
        #print Header+"** Warning: 3-ph fitting of the data failed for tau="+str(pulse_value*1E15)+" fs."
    
  #}}}
    
  print Header+"Generating comparison with Keldysh-Gruzdev model..." ## Looking forward to plot w_PI[intensity] from Keldysh_Gruzdev theory. 
  order = 50
  Efieldmin = np.min(Efield1_uniq)
  Efieldmax = np.max(Efield1_uniq)
  Efield1_Keldysh_log = np.linspace(np.log10(Efieldmin), np.log10(Efieldmax), 200)
  Efield1_Keldysh     = np.power(10., Efield1_Keldysh_log)
  print Efield1_Keldysh
  PeakWpi_Gruzdev     = GenerateKeldyshDatabase(Egap_SI, meff, pulseWavelength, Efield1_Keldysh, order)
  
  pplt.xlabel(r'Peak intensity (W/m$^2$)')
  pplt.ylabel(r'Final exc. density (m$^{-3}$)')
#  pplt.ylim((1E16,4.*5e28))
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  ax1.grid()
  pplt.tight_layout()
  if(PlotWpiAxis==True): #{{{
    ax12 = ax1.twinx()
    for pulse_value in pulses_uniq: #{{{
        # select the corresponding wavelength in results
        PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses, pulse_value, IndexTable['pulses1'])
        PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses_pulses, pulse_value, IndexTable['pulses2']) #keeps only pulses of the same duration.

        runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_pulses)
        
        lists = sorted(itertools.izip(*[FieldToIntensity(Efield1), Ne_final, pulses1]))
        Intensity_sorted, Ne_final_sorted, pulses1_sorted = list(itertools.izip(*lists))
        
        plot121, = ax12.loglog(Intensity_sorted, np.divide(Ne_final_sorted, np.multiply(2.,pulses1_sorted)), 'o--', label=r'$\tilde{w}_{PI}$, $\tau=$'+str(int(pulse_value*1E15))+' fs')
        plot1 += [plot121]
    #}}}
    plot122, = ax12.loglog(FieldToIntensity(Efield1_Keldysh), PeakWpi_Gruzdev, 'k-', label=r'$w_{PI}$ [Gruzdev], $\lambda=$'+str(round(pulseWavelength*1E9))+' nm')
    plot1    += [plot122]
    #ax12.set_yscale('log')
    #ax12.legend(loc='best', fancybox=True, framealpha=0)
    pplt.ylabel('$w_{PI}$ (m$^{-3}$.s$^{-1}$)')
    pplt.tight_layout()
  #}}}
    
  labels1 = [l.get_label() for l in plot1]
  ax1.legend(plot1, labels1, loc='best', fancybox=True, framealpha=0)
  #ax1.legend(loc=2, borderaxespad=0.) #bbox_to_anchor=(1.05, 1), 
  
  if(PlotVirtualExcitation):
    
    ax2 = pplt.subplot(plotB)
    plot1 = [] #np.empty((0,0))
    for pulse_value in pulses_uniq: #{{{  
        #pplt.title('Virtual exc. factor')
        # select the corresponding wavelength in results
        PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses, pulse_value, IndexTable['pulses1'])
        PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses_pulses, pulse_value, IndexTable['pulses2']) #keeps only pulses of the same duration.

        runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_pulses)
        
        print Header+"** Info: saving the data in CSV file for fitting. "
        filename = "EffectOfPulseIntensity-Duration-"+str(pulseWavelength*1E9)+"nm-tau"+str(int(pulse_value*1E15))+"fs"
        print Header+"** Info: saving in "+filename+".csv"
        #ExportToTxt(PlottingSinglePulses_pulses, filename+".csv")
        lists = sorted(itertools.izip(*[FieldToIntensity(Efield1), irreversibility, pulses1]))
        Intensity_sorted, irreversibility_sorted, pulses1_sorted = list(itertools.izip(*lists))
        
        plot_field, = ax2.plot(Intensity_sorted, irreversibility_sorted, '-o', label=r'$\zeta$, $\tau_1=$'+str(round(pulse_value*1E15))+'fs')
        plot1 += [plot_field]
    #}}}
    pplt.xlabel(r'Peak intensity (W/m$^2$)')
    pplt.ylabel(r'Virtual exc. $\zeta$')
    pplt.ylim((0,1))
    ax2.set_xscale('log')
    ax2.grid()
    
    ax22 = ax2.twinx()
    for pulse_value in pulses_uniq:
      # select the corresponding wavelength in results
      PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses, pulse_value, IndexTable['pulses1']) #keeps only 800 nm wavelenghts.

      runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_pulses)
      
      plot_field, = ax22.plot(FieldToIntensity(Efield1), gammaKeldysh(Egap_SI*e, meff, Efield1, wavelengths1), 'k-', label=r'$\gamma$')
      plot1 += [plot_field]
    pplt.ylabel(r'Adiabadicity $\gamma$')
    #pplt.ylim((0,1))
    ax22.set_yscale('log')
      
    # ax2.set_yscale("log")
    
    labels1 = [l.get_label() for l in plot1]
    ax22.legend(plot1, labels1, loc='best', fancybox=True, framealpha=0)
    #ax1.legend(loc=1, borderaxespad=0.)
    # ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #ax22.legend(loc=2, borderaxespad=0.)
    pplt.tight_layout()
  if(PlotEnergy):
    ax3 = pplt.subplot(plotC)
    #pplt.title('Virtual exc. factor')
    
    plot3 = []
    for pulse_value in pulses_uniq: #{{{
        # select the corresponding wavelength in results
        PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses, pulse_value, IndexTable['pulses1'])
        PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses_pulses, pulse_value, IndexTable['pulses2']) #keeps only pulses of the same duration.

        runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_pulses)
        
        lists = sorted(itertools.izip(*[FieldToIntensity(Efield1), Energy_max, pulses1]))
        Intensity_sorted, Energy_max_sorted, pulses1_sorted = list(itertools.izip(*lists))
        #print Header+"Ready to plot: "
        
        plot_field, = ax3.plot(Intensity_sorted, Energy_max_sorted-np.min(Energy_max_sorted)+EnergyMin, '-o', label=r'$\tau=$'+str(int(pulse_value*1E15))+" fs")
        plot3 += [plot_field]
    #}}}
    labels1 = [l.get_label() for l in plot3]
    ax3.legend(plot3, labels1, loc='best')
    pplt.xlabel(r'Peak intensity (W/m$^2$)')
    pplt.ylabel('Energy (eV)')
    #pplt.ylim((0,1))
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.grid()
    #ax32 = ax3.twinx()
    #for pulse_value in pulses_uniq:
      # select the corresponding wavelength in results
      #PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses, pulse_value, IndexTable['pulses1']) #keeps only 800 nm wavelenghts.

      #runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_pulses)
      
      #pplt.plot(FieldToIntensity(Efield1), gammaKeldysh(2.56*e, 0.228, Efield1, wavelengths1), '--', label=r'Adiabaticity, $\tau=$'+str(round(pulse_value*1E15))+' fs')
    #pplt.ylabel('Adiabadicity factor')
    #pplt.ylim((0,1))
    #ax32.set_yscale('log')
      
    # ax2.set_yscale("log")
    #ax3.legend(loc=1, borderaxespad=0.)
    # ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #ax32.legend(loc=2, borderaxespad=0.)
  
  
  # Saving the plot
  pplt.tight_layout()
  ExportName = "EffectOfPulseIntensity-Duration-"+str(pulseWavelength*1E9)+"nm"
  pplt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.eps')
  pplt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.png')
  print Header+"** Saved: 'Metaplot-SinglePulse-SingleColor-'"+ExportName+".png."
  if(ShowFinalPlot==1):
    pplt.show()

  del PlottingSinglePulses, wavelength_uniq, CEP1_uniq
  if(PlotVirtualExcitation): 
    del ax2, ax22
  if(PlotWpiAxis): 
    del ax12
  if(PlotEnergy): 
    del ax3
#}}}

## Single-pulse: Effects of CEP on excitation (for now), later: current and energy, and harmonic spectrum
def plotEffectOfIntensity_CEP(PlottingTable, pulseWavelength=800e-9, pulseDuration=20e-15, FilterCEP=False, CEP0=0., PulseShape=3, pulseRiseTime=0, ShowFinalPlot=0, Egap_SI=2.56*e):
  print "#====== Plot: SINGLE PULSE, EFFECT OF THE CEP, SAME WAVELENGTHS, SAME POLARIZATIONS: final density as function of field ======"
  #Filter: double pulses, same wavelengths, same polarization

  # Simple pulse, single color
  PlottingSinglePulses = PlottingTable #get the full database
  PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
  PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, pulseWavelength, IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0, IndexTable['FreezeHXC']) #keeps only non-frozen XC. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, PulseShape, IndexTable['PulseType']) #keeps only non-frozen XC. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, pulseDuration, IndexTable['pulses1']) #removes fields where tau != pulseDuration
  if(FilterCEP==True):
      PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, CEP0, IndexTable['CEP1']) #removes fields where tau != pulseDuration
  if(pulseRiseTime != 0): 
    PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, pulseRiseTime, IndexTable['RiseTime']) #removes fields where tau != pulseDuration

  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = CountTheNumberOfLaserParameters(PlottingSinglePulses)
  
  pplt.figure(figsize=(SizeX,SizeY))
  ax1 = pplt.subplot(111) #111 or plotA if using plotB
  pplt.title(r'$n_{exc}(I)$, $\tau=$'+str(int(pulseDuration*1E15))+' fs.')
  # Number of sweeped parameters: 3 (Efield, Delay, CrossPolarization)
  #pplt.scatter(FieldToIntensity(Efield1), Ne_final, c=Delay, s=(30.+CrossPolarization*50.), label=r'$N_{exc}$[Field]')
  plot1 = [] #np.empty((0,0))
  for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
    PlottingSinglePulses_CEP = FilterDatabase(PlottingSinglePulses, CEP1_value, IndexTable['CEP1']) #select the local CEP
    
    ## prints the current state of selection and stores the values of sweeped parameters
    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_CEP)
    
    #plot_field, = ax1.plot(FieldToIntensity(Efield1), Ne_final, '-', label=r'$N_{exc}^{final}$, $\phi_{CE}$='+str(round(CEP1_value*180./pi)))
    #lists = sorted(itertools.izip(*[FieldToIntensity(Efield1), Ne_final]))
    lists = sorted(itertools.izip(*[FieldToIntensity(Efield1), Ne_final, pulses1]))
    Intensity_sorted, Ne_final_sorted, pulses1_sorted = list(itertools.izip(*lists))
    #print Header+"Ready to plot: "
    #print "** Info: Sorting data..."
    #print Efield1_sorted, Ne_final_sorted
    plot_field, = ax1.plot(Intensity_sorted, Ne_final_sorted, 'o', label=r'$N_{exc}^{final}$, $\phi_{CE}$='+str(round(CEP1_value*180./pi)))
    print Header+"plotEffectOfIntensity_CEP(): adding CEP "+str(round(CEP1_value*180./pi))+" to the the plot."
    print plot_field
    print Header+"{Fore.RED}** Debug info: {Style.RESET_ALL}"+str(np.shape(plot_field))
    plot1 += [plot_field] # does cumulate each CEP now. Type: Line2D, not np.array. 
    
    #ax1.legend(loc=4, borderaxespad=0.) #bbox_to_anchor=(1.05, 1), 
  labels1 = [l.get_label() for l in plot1]
  ax1.legend(plot1, labels1, loc=2)
  plt.xlabel(r'Peak intensity (W/m$^2$)')
  plt.ylabel(r'Final density (m$^{-3}$)')
  #  pplt.ylim((1E16,4.*5e28))
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  #plot1_list = plot1.tolist()
    #print Efield1, CEP1_value
    
  #Adding axis of gammaKeldysh, but for all CEPs !!
  #runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses)
  
  ax12 = ax1.twinx()
  for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
    PlottingSinglePulses_CEP = FilterDatabase(PlottingSinglePulses, CEP1_value, IndexTable['CEP1']) #select the local CEP
    
    ## prints the current state of selection and stores the values of sweeped parameters
    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_CEP)
    
    #plot_field, = ax1.plot(FieldToIntensity(Efield1), Ne_final, '-', label=r'$N_{exc}^{final}$, $\phi_{CE}$='+str(round(CEP1_value*180./pi)))
    #lists = sorted(itertools.izip(*[FieldToIntensity(Efield1), Ne_final]))
    lists = sorted(itertools.izip(*[FieldToIntensity(Efield1), Ne_final, pulses1]))
    Intensity_sorted, Ne_final_sorted, pulses1_sorted = list(itertools.izip(*lists))
    plt.plot(Intensity_sorted, Ne_final_sorted/(pulses1), '-') #, label=r'$N_{exc}/\tau$') #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
  
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
  #ax1.legend(bbox_to_anchor=(1.15, 1), loc=2, borderaxespad=0.)
  #ax22.legend(loc=1, borderaxespad=0.)

  # Saving the plot
  ExportName = "EffectOfIntensityPerCEP"
  plt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.eps')
  plt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.png')
  if(ShowFinalPlot==1):
    pplt.show()

  del ax1, PlottingSinglePulses, wavelength_uniq
#}}}
  
def plotEffectOfCEP(PlottingTable, pulseWavelength=800e-9, pulseDuration=20e-15, FieldAmplitude=1E9, PulseShape=3, pulseRiseTime=0, ShowFinalPlot=0, Egap_SI=2.56*e):
  print "#{Fore.GREEN}====== Plot: SINGLE PULSE, EFFECT OF THE CEP, SAME WAVELENGTHS, SAME POLARIZATIONS: final density as function of field ======{Style.RESET_ALL}"
  #Filter: double pulses, same wavelengths, same polarization

  # Simple pulse, single color
  PlottingSinglePulses = PlottingTable #get the full database
  PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
  PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, pulseWavelength, IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts. #TODO: this field should be making a loop in available wavelengths given in wavelengths1_uniq.
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0, IndexTable['FreezeHXC']) #keeps only non-frozen XC. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, PulseShape, IndexTable['PulseType']) #keeps only non-frozen XC. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, pulseDuration, IndexTable['pulses1']) #removes fields where tau != pulseDuration
  PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, FieldAmplitude*1.01, IndexTable['Efield1'])
  PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, FieldAmplitude*0.99, IndexTable['Efield1'])
  if(pulseRiseTime != 0): 
    PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, pulseRiseTime, IndexTable['RiseTime']) #removes fields where tau != pulseDuration

  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = CountTheNumberOfLaserParameters(PlottingSinglePulses)
  
  pplt.figure(figsize=(SizeX,SizeY))
  ax1 = pplt.subplot(111) #111 or plotA if using plotB
  pplt.title(r'$n_{exc}(I)$, $\tau=$'+str(int(pulseDuration*1E15))+r' fs, $E_{peak}=$'+str(round(FieldAmplitude*1E-9))+' V/nm.')
  # Number of sweeped parameters: 3 (Efield, Delay, CrossPolarization)
  #pplt.scatter(FieldToIntensity(Efield1), Ne_final, c=Delay, s=(30.+CrossPolarization*50.), label=r'$N_{exc}$[Field]')
  plot1 = [] #np.empty((0,0))
  #for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
  #PlottingSinglePulses_CEP = FilterDatabase(PlottingSinglePulses, CEP1_value, IndexTable['CEP1']) #select the local CEP
  
  ## prints the current state of selection and stores the values of sweeped parameters
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses)
  
  #plot_field, = ax1.plot(FieldToIntensity(Efield1), Ne_final, '-', label=r'$N_{exc}^{final}$, $\phi_{CE}$='+str(round(CEP1_value*180./pi)))
  #lists = sorted(itertools.izip(*[FieldToIntensity(Efield1), Ne_final]))
  lists = sorted(itertools.izip(*[CEP1, Ne_final, pulses1]))
  CEP1_sorted, Ne_final_sorted, pulses1_sorted = list(itertools.izip(*lists))
  #print Header+"Ready to plot: "
  #print "** Info: Sorting data..."
  #print Efield1_sorted, Ne_final_sorted
  plot_field, = ax1.plot(np.multiply(CEP1_sorted,180./pi), Ne_final_sorted, 'o', label=r'$N_{exc}^{final}$')
  #print Header+"plotEffectOfIntensity_CEP(): adding CEP "+str(round(CEP1_value*180./pi))+" to the the plot."
  print plot_field
  print Header+"{Fore.RED}** Debug info: {Style.RESET_ALL}"+str(np.shape(plot_field))
  plot1 += [plot_field] # does cumulate each CEP now. Type: Line2D, not np.array. 
    
    #ax1.legend(loc=4, borderaxespad=0.) #bbox_to_anchor=(1.05, 1), 
  labels1 = [l.get_label() for l in plot1]
  ax1.legend(plot1, labels1, loc='best', fancybox=True, framealpha=0)
  #plt.xlabel(r'Peak intensity (W/m$^2$)')
  plt.xlabel(r'Carrier envelope phase $\phi$ (deg)')
  plt.ylabel(r'Final density (m$^{-3}$)')
  #  pplt.ylim((1E16,4.*5e28))
  #ax1.set_xscale('log')
  #ax1.set_yscale('log')
  #plot1_list = plot1.tolist()
    #print Efield1, CEP1_value
    
  #Adding axis of gammaKeldysh, but for all CEPs !!
  #runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses)
  
  ax12 = ax1.twinx()
  #for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
  #PlottingSinglePulses_CEP = FilterDatabase(PlottingSinglePulses, CEP1_value, IndexTable['CEP1']) #select the local CEP
  
  ## prints the current state of selection and stores the values of sweeped parameters
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses)
  
  #plot_field, = ax1.plot(FieldToIntensity(Efield1), Ne_final, '-', label=r'$N_{exc}^{final}$, $\phi_{CE}$='+str(round(CEP1_value*180./pi)))
  #lists = sorted(itertools.izip(*[FieldToIntensity(Efield1), Ne_final]))
  lists = sorted(itertools.izip(*[CEP1, Ne_final, pulses1]))
  CEP1_sorted, Ne_final_sorted, pulses1_sorted = list(itertools.izip(*lists))
  plt.plot(np.multiply(CEP1_sorted,180./pi), np.divide(Ne_final_sorted, pulses1_sorted), '-') #, label=r'$N_{exc}/\tau$') #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
  
  #ax12.set_yscale('log')
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
  #ax1.legend(bbox_to_anchor=(1.15, 1), loc=2, borderaxespad=0.)
  #ax22.legend(loc=1, borderaxespad=0.)

  # Saving the plot
  plt.tight_layout()
  ExportName = "EffectOfCEP"
  plt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.eps')
  plt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.png')
  if(ShowFinalPlot==1):
    pplt.show()

  del ax1, PlottingSinglePulses, wavelength_uniq
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
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_CEP_NoFreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), Ne_final, 'ow', label=r'$N_{exc}^{final}$') #, CEP='+str(round(CEP1_value*180./pi)))
  
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_CEP_FreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), Ne_final, 'ok', label=r'$N_{exc}^{final}$, FrozenHXC')#, CEP='+str(round(CEP1_value*180./pi)))
  
  ax1.legend(loc='best', fancybox=True, framealpha=0) #bbox_to_anchor=(1.05, 1), 
  #pplt.xlabel(r'Peak intensity (W/m$^2$)')
  pplt.ylabel(r'Final density (m$^{-3}$)')
  #  pplt.ylim((1E16,4.*5e28))
  ax1.set_xscale('log')
  ax1.set_yscale('log')

  #Adding axis of gammaKeldysh
  ax12 = ax1.twinx()
  
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_CEP_NoFreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), Ne_final/pulses1, 'ow') #, label=r'$N_{exc}/\tau$') #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
  
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_CEP_FreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), Ne_final/pulses1, 'ok') #, label=r'$N_{exc}/\tau$, FrozenHXC') #, I='+str(FieldToIntensity(Efield1_value))+r'
  
  ax12.set_yscale('log')
  ax12.set_ylabel(r'$w_{PI}$ (m$^{-3}$/s)')

  ax2 = pplt.subplot(plotB)
  pplt.title(r'Virtual excitation $\zeta$')
  
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_CEP_NoFreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), irreversibility, 'ow', label='Virtual exc.')
  
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_CEP_FreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), irreversibility, 'ok', label='Virtual exc., FreezeHXC')
  
  pplt.xlabel(r'Peak intensity (W/m$^2$)')
  pplt.ylabel(r'Virtual exc. $\zeta$')
  pplt.ylim((0,1))
  ax2.set_xscale('log')
  #adding the Keldysh adiabatic factor
  ax22 = ax2.twinx()
  
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_CEP_NoFreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), gammaKeldysh(Egap_SI*e, meff, Efield1, wavelengths1), 'or', label=r'Adiabaticity $\gamma$') #BUG: 2.56 must not be specified by hand here! 
  
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_CEP_FreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), gammaKeldysh(Egap_SI*e, meff, Efield1, wavelengths1), 'or') #, label='Adiabaticity') #BUG: 2.56 must not be specified by hand here! 
  
  pplt.ylabel(r'Adiabadicity $\gamma$')
  #pplt.ylim((0,1))
  ax22.set_yscale('log')

  #print Efield1, CEP1_value
  # ax2.set_yscale("log")
  #ax2.legend(loc=2, borderaxespad=0.)
  # ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
  #ax22.legend(loc=1, borderaxespad=0.)
  pplt.tight_layout()
  # Saving the plot
  pplt.savefig('Metaplot-SinglePulse-SingleColor-EffectOfFrozenHXC.eps')
  pplt.savefig('Metaplot-SinglePulse-SingleColor-EffectOfFrozenHXC.png')
  #if(ShowFinalPlot==1):
    #pplt.show()

  del ax1, PlottingSinglePulses, wavelength_uniq, wavelengths1
  #}}}

#del PlottingTable

## Plots the effect of wavelength (1D plot) by introducing 10 percents tolerance on the electric field
def plotEffectOfWavelength_Tolerance(PlottingTable, FieldTolerance=0.05, pulseDuration=20e-15, CEP=0., PulseType=3, RiseTime=2, ShowFinalPlot=0, Egap_SI=2.56*e):
  
  meff = 0.2226 #effective mass of electron
  
  print "#====== Plot: SINGLE PULSE, EFFECT OF THE WAVELENGTH, SAME CEP, SAME POLARIZATIONS: final density as function of field ======"
  #Filter: double pulses, same wavelengths, same polarization
  
  
  
  # Simple pulse, single color
  PlottingSinglePulses = PlottingTable #get the full database
  PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
  PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['CrossPolarization']) #removes fields where CrossPolarization is different from 0. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, PulseType, IndexTable['PulseType']) #keeps only pulses of type StH (3)
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, RiseTime, IndexTable['RiseTime']) #keeps only pulses of type STH RiseTime
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, pulseDuration, IndexTable['pulses1']) #keeps only pulses of selected duration in argument
  if(RiseTime >= 1): #select by RiseTime only if it is defined
      PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, RiseTime, IndexTable['RiseTime'])
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, CEP, IndexTable['CEP1']) #keeps only CEP = 0. 
  PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, 2100e-9, IndexTable['wavelengths1'])

  ## TODO: check if CEP could count at smaller wavelength. 

  print Header+"Available data: "
  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = CountTheNumberOfLaserParameters(PlottingSinglePulses)
  #if(len(CEP1_uniq) == 1): #{{{ #single color simulations
  # Complexity level: 2
  
  pplt.figure(figsize=(SizeX,SizeY)) #new figure
  plotA=211; plotB=212
  ax1 = pplt.subplot(plotA)
  #pplt.title('Final density as function of laser wavelength')
  print Header+"Total number of datasets: "+str(len(set(Efield1_uniq)))
  
  plot1 = [] #np.empty((0,0))
  for Efield1_value in Efield1_uniq: #for each value of Efield1, we make a new search, and plot one new line in the chart
    #PlottingSinglePulses_Efield = FilterDatabase(PlottingSinglePulses, Efield1_value, IndexTable['Efield1']) #select the local Efield
    # TODO: how to merge data with values close to each others? Maybe introduce tolerance of 10% ?
    PlottingSinglePulses_Efield = FilterDatabaseLowerThan(PlottingSinglePulses, Efield1_value*(1.+FieldTolerance), IndexTable['Efield1'])
    PlottingSinglePulses_Efield = FilterDatabaseGreaterThan(PlottingSinglePulses_Efield, Efield1_value*(1.-FieldTolerance), IndexTable['Efield1'])
    
    #print PlottingSinglePulses_Efield
    ## prints the current state of selection and stores the values of sweeped parameters
    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_Efield)
    
    lists = sorted(itertools.izip(*[wavelengths1, Ne_final, Energy_max])) #NOTE: Ne_max makes more sense than Ne_final
    wavelengths1_sorted, Ne_final_sorted, Energy_max_sorted = list(itertools.izip(*lists))
    wavelengths1_sorted = np.multiply(1e9, wavelengths1_sorted)
    #print Efield1_value
    # TODO: PUT here a select case which attributes a new color to a new field interval
    if(1e8*(1.- FieldTolerance) < Efield1_value < 1e8*(1.+ FieldTolerance)):
      #print Header+"Dataset for the field: "+str(Efield1_value*1e-9)+' V/nm'
      plot_field, = ax1.plot(wavelengths1_sorted, Ne_final_sorted, '^r', label=r'$E_{peak}=$'+str(round(Efield1_value*1E-9))+' V/nm')
      plot1 += [plot_field]
      #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
    if(1e9*(1.- FieldTolerance) < Efield1_value < 1e9*(1.+ FieldTolerance)):
      #print Header+"Dataset for the field: "+str(Efield1_value*1e-9)+' V/nm'
      plot_field, = ax1.plot(wavelengths1_sorted, Ne_final_sorted, '^b', label=r'$E_{peak}=$'+str(round(Efield1_value*1E-9))+' V/nm') 
      plot1 += [plot_field]
      #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
    if(5e9*(1.- FieldTolerance) < Efield1_value < 5e9*(1.+ FieldTolerance)):
      #print Header+"Dataset for the field: "+str(Efield1_value*1e-9)+' V/nm'
      plot_field, = ax1.plot(wavelengths1_sorted, Ne_final_sorted, '^g', label=r'$E_{peak}=$'+str(round(Efield1_value*1E-9))+' V/nm')
      plot1 += [plot_field]
      #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
    if(1e10*(1.- FieldTolerance) < Efield1_value < 1e10*(1.+ FieldTolerance)):
      #print Header+"Dataset for the field: "+str(Efield1_value*1e-9)+' V/nm'
      plot_field, = ax1.plot(wavelengths1_sorted, Ne_final_sorted, '^k', label=r'$E_{peak}=$'+str(round(Efield1_value*1E-9))+' V/nm')
      plot1 += [plot_field]
      #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
    del PlottingSinglePulses_Efield
    
  labels1 = [l.get_label() for l in plot1]
  #ax1.legend(plot1, labels1, loc='best')    
  #ax1.legend(loc=3, borderaxespad=0., bbox_to_anchor=(1.05, 0.05))
  #pplt.xlabel('Laser wavelength (nm)')
  pplt.ylabel(r'Final density $n_{exc}$ (m$^{-3}$)')
  pplt.xlim((200.,2100.))
  #ax1.set_xscale('log')
  ax1.set_yscale('log')
  ax1.grid()
  
  
  ax12 = ax1.twinx()
  for Efield1_value in Efield1_uniq: #for each value of Efield1, we make a new search, and plot one new line in the chart
    #PlottingSinglePulses_Efield = FilterDatabase(PlottingSinglePulses, Efield1_value, IndexTable['Efield1']) #select the local Efield
    PlottingSinglePulses_Efield = FilterDatabaseLowerThan(PlottingSinglePulses, Efield1_value*(1.+FieldTolerance), IndexTable['Efield1'])
    PlottingSinglePulses_Efield = FilterDatabaseGreaterThan(PlottingSinglePulses_Efield, Efield1_value*(1.-FieldTolerance), IndexTable['Efield1'])
    #print Header+"Dataset for the field: "+str(Efield1_value*1e-9)+' V/nm'
    #print PlottingSinglePulses_Efield
    # TODO: as fields are all slighly different, we should project on the neareast neighboring value, up to 1E8 V/m. 
    
    ## prints the current state of selection and stores the values of sweeped parameters
    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_Efield)
    
    lists = sorted(itertools.izip(*[wavelengths1, Ne_final, pulses1, Energy_max])) #NOTE: Ne_max makes more sense than Ne_final
    wavelengths1_sorted, Ne_final_sorted, pulses1_sorted, Energy_max_sorted = list(itertools.izip(*lists))
    
    wavelengths1_sorted_show = np.multiply(1e9, wavelengths1_sorted)
    wpi_sorted = np.divide(Ne_final_sorted, np.multiply(2e0, pulses1_sorted))
    
    print Header+"Building comparison with Keldysh-Gruzdev..."
    order = 50 
    wavelengths1_min = np.min(wavelengths1); wavelengths1_max = np.max(wavelengths1)
    wavelengths1_Keldysh = np.linspace(wavelengths1_min, wavelengths1_max, 100)
    wavelengths1_Keldysh_show = np.multiply(1E9, wavelengths1_Keldysh)
    PeakWpi_Gruzdev = GenerateKeldyshDatabase(Egap_SI, meff, wavelengths1_Keldysh, Efield1_value, order)
    print PeakWpi_Gruzdev
    
    #print Efield1_value
    # TODO: PUT here a select case which attributes a new color to a new field interval
    if(1e8*(1.-FieldTolerance) < Efield1_value < 1e8*(1.+FieldTolerance)):
      ax12.plot(wavelengths1_sorted_show, wpi_sorted,       '+-r', label=r'$\widetilde{w}_{PI}$, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
      ax12.plot(wavelengths1_Keldysh_show, PeakWpi_Gruzdev, '--r', label=r'$\widetilde{w}_{PI}$ [Gruzdev], I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
    if(1e9*(1.-FieldTolerance) < Efield1_value < 1e9*(1.+FieldTolerance)):
      ax12.plot(wavelengths1_sorted_show, wpi_sorted,       '+-b', label=r'$\widetilde{w}_{PI}$, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
      ax12.plot(wavelengths1_Keldysh_show, PeakWpi_Gruzdev, '--b', label=r'$\widetilde{w}_{PI}$ [Gruzdev], I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
    if(5e9*(1.-FieldTolerance) < Efield1_value < 5e9*(1.+FieldTolerance)):
      ax12.plot(wavelengths1_sorted_show, wpi_sorted,       '+-g', label=r'$\widetilde{w}_{PI}$') #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
      ax12.plot(wavelengths1_Keldysh_show, PeakWpi_Gruzdev, '--g', label=r'$\widetilde{w}_{PI}$ [Gruzdev], I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
    if(1e10*(1.-FieldTolerance) < Efield1_value < 1e10*(1.+FieldTolerance)):
      ax12.plot(wavelengths1_sorted_show, wpi_sorted,       '+-k', label=r'$\widetilde{w}_{PI}$') #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
      ax12.plot(wavelengths1_Keldysh_show, PeakWpi_Gruzdev, '--k', label=r'$\widetilde{w}_{PI}$ [Gruzdev], I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
    del PlottingSinglePulses_Efield
  ax12.set_yscale('log')
  ax12.set_ylabel(r'$w_{PI}$ (m$^{-3}$/s)')

  ax2 = plt.subplot(plotB)
  #plt.title('Irreversibility factor')
  #pplt.plot(FieldToIntensity(Efield1), irreversibility, 'o', label='Irreversibility')
  for Efield1_value in Efield1_uniq: #for each value of Efield1, we make a new search, and plot one new line in the chart
    #PlottingSinglePulses_Efield = FilterDatabase(PlottingSinglePulses, Efield1_value, IndexTable['Efield1']) #select the local Efield
    PlottingSinglePulses_Efield = FilterDatabaseLowerThan(PlottingSinglePulses, Efield1_value*(1.+FieldTolerance), IndexTable['Efield1'])
    PlottingSinglePulses_Efield = FilterDatabaseGreaterThan(PlottingSinglePulses_Efield, Efield1_value*(1.-FieldTolerance), IndexTable['Efield1'])
    #print Header+"Dataset for the field: "+str(Efield1_value*1e-9)+' V/nm'
    #print PlottingSinglePulses_Efield
    # TODO: as fields are all slighly different, we should project on the neareast neighboring value, up to 1E8 V/m. 
    
    ## prints the current state of selection and stores the values of sweeped parameters
    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_Efield)
    
    lists = sorted(itertools.izip(*[wavelengths1, Ne_final, Energy_max])) #NOTE: Ne_max makes more sense than Ne_final
    wavelengths1_sorted, Ne_final_sorted, Energy_max_sorted = list(itertools.izip(*lists))
    wavelengths1_sorted = np.multiply(1e9, wavelengths1_sorted)
    EnergyMin = 1e-5 #0e0 #1E-7 #shifting energies to avoid zero on logscale. 
    #print Efield1_value
    #TODO: PUT here a select case which attributes a new color to a new field interval
    if(1e8*(1.-FieldTolerance) < Efield1_value < 1e8*(1.+FieldTolerance)):
      ax2.plot(wavelengths1_sorted, Energy_max_sorted-np.min(Energy_max_sorted)+EnergyMin, '-or')
    if(1e9*(1.-FieldTolerance) < Efield1_value < 1e9*(1.+FieldTolerance)):
      ax2.plot(wavelengths1_sorted, Energy_max_sorted-np.min(Energy_max_sorted)+EnergyMin, '-ob')
    if(5e9*(1.-FieldTolerance) < Efield1_value < 5e9*(1.+FieldTolerance)):
      ax2.plot(wavelengths1_sorted, Energy_max_sorted-np.min(Energy_max_sorted)+EnergyMin, '-og')
    if(1e10*(1.-FieldTolerance) < Efield1_value < 1e10*(1.+FieldTolerance)):
      ax2.plot(wavelengths1_sorted, Energy_max_sorted-np.min(Energy_max_sorted)+EnergyMin, '-ok')
    del PlottingSinglePulses_Efield
    
  #pplt.xlabel(r'Intensity (W/m$^2$)')
  pplt.xlabel('Laser wavelength (nm)')
  #pplt.ylabel('Irreversibility factor')
  pplt.ylabel('Energy (eV)')
  #pplt.ylim((0,1))
  #ax2.set_xscale('log')
  ax22 = ax2.twinx()
  
  for Efield1_value in Efield1_uniq: #for each value of Efield1, we make a new search, and plot one new line in the chart
    #PlottingSinglePulses_Efield = FilterDatabase(PlottingSinglePulses, Efield1_value, IndexTable['Efield1']) #select the local Efield
    PlottingSinglePulses_Efield = FilterDatabaseLowerThan(PlottingSinglePulses, Efield1_value*(1.+FieldTolerance), IndexTable['Efield1'])
    PlottingSinglePulses_Efield = FilterDatabaseGreaterThan(PlottingSinglePulses_Efield, Efield1_value*(1.-FieldTolerance), IndexTable['Efield1'])
    #print Header+"Dataset for the field: "+str(Efield1_value*1e-9)+' V/nm'
    #print PlottingSinglePulses_Efield
    # TODO: as fields are all slighly different, we should project on the neareast neighboring value, up to 1E8 V/m. 
    
    ## prints the current state of selection and stores the values of sweeped parameters
    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_Efield)
    
    lists = sorted(itertools.izip(*[wavelengths1, Ne_final, Energy_max, irreversibility])) #NOTE: Ne_max makes more sense than Ne_final
    wavelengths1_sorted, Ne_final_sorted, Energy_max_sorted, irreversibility_sorted = list(itertools.izip(*lists))
    wavelengths1_sorted_show = np.multiply(1e9, wavelengths1_sorted)
    
    #print Efield1_value
    # TODO: PUT here a select case which attributes a new color to a new field interval
    if(1e8*(1.-FieldTolerance) < Efield1_value < 1e8*(1.+FieldTolerance)):
      ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+r', label='Irreversibility')
      #ax22.plot(wavelengths1_sorted_show, gammaKeldysh(2.56*e, 0.228, Efield1_value, wavelengths1_sorted), '--r', label='Adiabaticity')
    if(1e9*(1.-FieldTolerance) < Efield1_value < 1e9*(1.+FieldTolerance)):
      ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+b', label='Irreversibility')
      #ax22.plot(wavelengths1_sorted_show, gammaKeldysh(2.56*e, 0.228, Efield1_value, wavelengths1_sorted), '--b', label='Adiabaticity')
    if(5e9*(1.-FieldTolerance) < Efield1_value < 5e9*(1.+FieldTolerance)):
      ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+g', label='Irreversibility')
      #ax22.plot(wavelengths1_sorted_show, gammaKeldysh(2.56*e, 0.228, Efield1_value, wavelengths1_sorted), '--g', label='Adiabaticity')
    if(1e10*(1.-FieldTolerance) < Efield1_value < 1e10*(1.+FieldTolerance)):
      ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+k', label='Irreversibility')
      #ax22.plot(wavelengths1_sorted_show, gammaKeldysh(2.56*e, 0.228, Efield1_value, wavelengths1_sorted), '--k', label='Adiabaticity')
    del PlottingSinglePulses_Efield
  
  
  pplt.ylabel(r'Irreversibility $\zeta$')
  #pplt.ylabel(r'Adiabaticity $\gamma$, Irreversibility $\zeta$')
  #ax22.set_ylim((0,1))
  ax2.set_yscale("log")
  ax22.set_yscale('log')
  #ax22.grid()
  pplt.xlim((200.,2100.))
  
  ResonantWavelength = h*c/Egap_SI
  print "** Adding vertical bar for the resonant transitions at: "+str(round(1e9*ResonantWavelength))+"nm."
  
  ax1.axvline(x=1e9*ResonantWavelength, color='k', linestyle='--', linewidth=1)
  ax1.axvline(x=1e9*2*ResonantWavelength, color='k', linestyle='--', linewidth=1)
  ax1.axvline(x=1e9*3*ResonantWavelength, color='k', linestyle='--', linewidth=1)
  ax1.axvline(x=1e9*4*ResonantWavelength, color='k', linestyle='--', linewidth=1)
  
  ax2.axvline(x=1e9*ResonantWavelength, color='k', linestyle='--', linewidth=1)
  ax2.axvline(x=1e9*2*ResonantWavelength, color='k', linestyle='--', linewidth=1)
  ax2.axvline(x=1e9*3*ResonantWavelength, color='k', linestyle='--', linewidth=1)
  ax2.axvline(x=1e9*4*ResonantWavelength, color='k', linestyle='--', linewidth=1)
  #ax2.legend(loc=2, borderaxespad=0.)
  #ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
  #ax22.legend(loc=1, borderaxespad=0.)

  # Saving the plot
  pplt.savefig('Metaplot-SinglePulse-SingleColor-EffectOfWavelength_Tol10.eps')
  pplt.savefig('Metaplot-SinglePulse-SingleColor-EffectOfWavelength_Tol10.png')
  if(ShowFinalPlot==1):
    pplt.show()

  del ax1, PlottingSinglePulses
  #}}}

## Plots the N_exc(wavelength, field) (2D plot)
# VALID, but of course interference of 2 waves is oscillating very much !
def plotEffectOfWavelength_2dMap(PlottingTable, pulseDuration=20e-15, PulseType=2, RiseCycles=0, ShowFinalPlot=0, Egap_SI=2.56*e):
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
      runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_Efield_wavelength)
      
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
  pplt.tight_layout()
  pplt.savefig('Metaplot-SinglePulse-SingleColor-EffectOfWavelengthAndField.eps')
  pplt.savefig('Metaplot-SinglePulse-SingleColor-EffectOfWavelengthAndField.png')
  if(ShowFinalPlot==1):
    pplt.show()

  del PlottingSinglePulses, DensityMax2D
#}}}


# Plots Ne_final(Efield1, Efield2) for pulses of fixed wavelengths, same pulse duration, same CEP
# THIS FIGURE IS THE ONE THAT CAN BE COMPARED WITH VLADIMIR ZHUKOV
def plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1=800e-9, wavelength2=1600e-9, PulseDuration=20e-15, PulseType=3, RiseCycles=2, WhatToPlot=0, Egap_SI=2.56*e): #{{{
  print "====== Plot: DOUBLE PULSE, SAME POLARIZATION, SAME PULSE DURATIONS, Ne(E1, E2) ====="
  #Filters: double pulses, same wavelengths, same polarization
  PlottingDoublePulses = PlottingTable #get the full database
  PlottingDoublePulses = FilterDatabaseGreaterThan(PlottingDoublePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak < 1 V/m. #works
  PlottingDoublePulses = FilterDatabaseGreaterThan(PlottingDoublePulses, 1E0, IndexTable['Efield2']) #removes fields where Epeak2 < 1 V/m to select double pulses.
  PlottingDoublePulses = FilterDatabase(PlottingDoublePulses, float(0.), IndexTable['CrossPolarization']) #removes fields where Polarization1 != Polarization2
  PlottingDoublePulses = FilterDatabase(PlottingDoublePulses, PulseDuration, IndexTable['pulses2']) # Removing fields were pulseDuration1 != pulseDuration2
  PlottingDoublePulses = FilterDatabase(PlottingDoublePulses, wavelength1, IndexTable['wavelengths1']) # Selecting color1
  PlottingDoublePulses = FilterDatabase(PlottingDoublePulses, wavelength2, IndexTable['wavelengths2']) # Selecting color2
  PlottingDoublePulses = FilterDatabase(PlottingDoublePulses, PulseType, IndexTable['PulseType']) # Selecting pulse type
  if(RiseCycles != 0): #0 selects all the possible types
    PlottingDoublePulses = FilterDatabase(PlottingDoublePulses, RiseCycles, IndexTable['RiseTime']) # Selecting Rise Time
  #### 2D PLOTTING Density(Efield1, Efield2) for one given set of laser parameters. For each combination, we produce one plot like this. 
  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq, wavelength2_uniq, CEP2_uniq, Efield2_uniq, pulses2_uniq= CountTheNumberOfLaserParameters(PlottingDoublePulses)
  print Header+"Total number of datasets: "+str(len(set(Efield1_uniq))*len(set(Efield2_uniq)))
  
  TotalField2D   = np.empty((0,3)) #preparing the 2D plot
  DensityMax2D = np.empty((0,3)) #preparing the 2D plot
  CurrentMax2D = np.empty((0,3)) #preparing the 2D plot
  EnergyMax2D  = np.empty((0,3)) #preparing the 2D plot
  WpiMax2D     = np.empty((0,3)) #preparing the 2D plot
  for Efield1_i in Efield1_uniq: #{{{ #selecting amplitude1
    PlottingDoublePulses_SelectedField1 = FilterDatabase(PlottingDoublePulses, Efield1_i, IndexTable['Efield1'])
    for Efield2_i in Efield2_uniq: #{{{ #selecting amplitude2
      # Selecting Field2
      PlottingDoublePulses_SelectedField2 = FilterDatabase(PlottingDoublePulses_SelectedField1, Efield2_i, IndexTable['Efield2'])
      
      # Constructing the density table
      runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingDoublePulses_SelectedField2)
      
      print Header+"Field1: "+str(Efield1)+" m, Field2: "+str(Efield2)+" m."
      if(len(Ne_final) >= 1):
        ToBeAddedField  =np.array([Efield1[0], Efield2[0], Efield1[0]+Efield2[0]]) #builds the total peak field
        ToBeAdded       =np.array([Efield1[0], Efield2[0], Ne_final[0]]) #BUG: removes some results, not clean enough!
        ToBeAddedCurrent=np.array([Efield1[0], Efield2[0], TotalCurrent_final[0]]) #BUG: removes some results, not clean enough!
        ToBeAddedEnergy =np.array([Efield1[0], Efield2[0], Energy_max[0]]) #BUG: removes some results, not clean enough!
        ToBeAddedWpi    =np.array([Efield1[0], Efield2[0], W_pi[0]]) #BUG: removes some results, not clean enough!
        print ToBeAdded
        TotalField2D=np.vstack((TotalField2D, ToBeAddedField))
        DensityMax2D=np.vstack((DensityMax2D, ToBeAdded))
        CurrentMax2D=np.vstack((CurrentMax2D, ToBeAddedCurrent))
        EnergyMax2D =np.vstack((EnergyMax2D, ToBeAddedEnergy))
        WpiMax2D    =np.vstack((WpiMax2D   , ToBeAddedWpi))
        print "Cleaning..."
        del ToBeAdded, ToBeAddedCurrent, ToBeAddedEnergy, ToBeAddedWpi, ToBeAddedField
    #}}}
  #}}}
  PulseTypeName = np.asfarray(PulseType)
  # Stack is built. Export, show and plot on a heat map. 
  print "** Info: exporting the sub-plot2D: Ne_final(Field1, Field2)..."
  print "** Info: useful to visualize available data with GNUplot. "
  filename_field   = "Bicolor-"+str(PulseTypeName)+"-Color1-"+str(round(wavelength1*1E9))+"-Color2-"+str(round(wavelength2*1E9))+"-"+str(int(PulseDuration*1e15))+"fs-EPeakWithE1E2"
  filename         = "Bicolor-"+str(PulseTypeName)+"-Color1-"+str(round(wavelength1*1E9))+"-Color2-"+str(round(wavelength2*1E9))+"-"+str(int(PulseDuration*1e15))+"fs-NeFinalWithE1E2"
  filename_current = "Bicolor-"+str(PulseTypeName)+"-Color1-"+str(round(wavelength1*1E9))+"-Color2-"+str(round(wavelength2*1E9))+"-"+str(int(PulseDuration*1e15))+"fs-CurrentFinalWithE1E2"
  filename_energy  = "Bicolor-"+str(PulseTypeName)+"-Color1-"+str(round(wavelength1*1E9))+"-Color2-"+str(round(wavelength2*1E9))+"-"+str(int(PulseDuration*1e15))+"fs-EnergyFinalWithE1E2"
  filename_wpi     = "Bicolor-"+str(PulseTypeName)+"-Color1-"+str(round(wavelength1*1E9))+"-Color2-"+str(round(wavelength2*1E9))+"-"+str(int(PulseDuration*1e15))+"fs-WpiWithE1E2"
  # Sorting the data before storing and plotting
  TotalField2Dsorted = TotalField2D[np.lexsort((TotalField2D[:,2], TotalField2D[:,1], TotalField2D[:,0]))]
  DensityMax2Dsorted = DensityMax2D[np.lexsort((DensityMax2D[:,2], DensityMax2D[:,1], DensityMax2D[:,0]))]
  CurrentMax2Dsorted = CurrentMax2D[np.lexsort((CurrentMax2D[:,2], CurrentMax2D[:,1], CurrentMax2D[:,0]))]
  EnergyMax2Dsorted  = EnergyMax2D[np.lexsort((EnergyMax2D[:,2], EnergyMax2D[:,1], EnergyMax2D[:,0]))]
  WpiMax2Dsorted     = WpiMax2D[np.lexsort((WpiMax2D[:,2], WpiMax2D[:,1], WpiMax2D[:,0]))]
  
  
  ## Lazy to build up a new function for that for now.
  if(WhatToPlot == 0): 
      QuantityTitle=r'$N_{exc}(E_1,E_2)$, $(\lambda_1, \lambda_2) = $('+str(wavelength1*1E9)+', '+str(wavelength2*1E9)+') nm'
  elif(WhatToPlot == 4): 
      DensityMax2Dsorted = TotalField2Dsorted; 
      filename = filename_field
      QuantityTitle=r"$E_{tot}(E_1,E_2)$, $(\lambda_1, \lambda_2) = $("+str(wavelength1*1E9)+", "+str(wavelength2*1E9)+") nm."
  elif(WhatToPlot == 1): 
      DensityMax2Dsorted = CurrentMax2Dsorted; 
      filename = filename_current
      QuantityTitle=r"$J_{end}(E_1,E_2)$, $(\lambda_1, \lambda_2) = $("+str(wavelength1*1E9)+", "+str(wavelength2*1E9)+") nm."
  elif(WhatToPlot == 2): 
      DensityMax2Dsorted = EnergyMax2Dsorted
      filename = filename_energy
      QuantityTitle=r"$\xi_{end}(E_1,E_2)$, $(\lambda_1, \lambda_2) = $("+str(wavelength1*1E9)+", "+str(wavelength2*1E9)+") nm."
  elif(WhatToPlot == 3): 
      DensityMax2Dsorted = WpiMax2Dsorted
      filename = filename_wpi
      QuantityTitle=r'$w_{PI}(E_1,E_2)$, $(\lambda_1, \lambda_2) = $('+str(wavelength1*1E9)+', '+str(wavelength2*1E9)+') nm'
  else: 
      print Header+"** plot_Zhukov_BicolorMaps: wrong request!"
  
  # Outputing the data to plot from outside
  ExportToTxt(TotalField2Dsorted, filename+".csv")
  ExportToTxt(DensityMax2Dsorted, filename+".csv")
  ExportToTxt(CurrentMax2Dsorted, filename+".csv")
  ExportToTxt(EnergyMax2Dsorted, filename+".csv")
  ExportToTxt(WpiMax2Dsorted, filename+".csv")
  
  # Shaping the tables
  BrutX = 1E-9*DensityMax2Dsorted[:,0]
  BrutY = 1E-9*DensityMax2Dsorted[:,1]
  BrutZ = DensityMax2Dsorted[:,2]
  
  print len(BrutX)
  print len(BrutY)
  print len(BrutZ)
  X=list(set(BrutX))
  Y=list(BrutY)
  Z=DensityMax2Dsorted
  YY, XX = np.meshgrid(X, Y)
  # If the size is big enough, let's make a plot.
  if(len(DensityMax2D) >= 2): #{{{
    try:
        triang = tri.Triangulation(BrutX, BrutY)
        triang.set_mask(np.hypot(BrutX[triang.triangles].mean(axis=1),
                                BrutY[triang.triangles].mean(axis=1))
                        < 0.15)
        refiner = tri.UniformTriRefiner(triang)
        tri_refi, z_test_refi = refiner.refine_field(BrutZ, subdiv=3)
        #plt.figure()
        #plt.title(r'$N_{exc}(E_1,E_2)$, $\lambda_1=$'+str(wavelength1*1E9)+r'nm, $\lambda_2=$'+str(wavelength2*1E9)+'nm.')
        fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')
        #plt.gca().set_aspect('equal')
        plt.title(QuantityTitle)
        plt.xlabel(r'$E_1$ (V/nm)')
        plt.ylabel(r'$E_2$ (V/nm)')
        #CS = plt.imshow(XX, YY, DensityMax2Dsorted) #, cmap=plt.cm.Blues)
        #CS = plt.contourf(XX, YY, DensityMax2Dsorted, cmap=plt.cm.Blues)
        #CS = ax.scatter(BrutX, BrutY, BrutZ) #, cmap=plt.cm.Blues) #WORKS but limited to boring dots
        #CS = ax.plot_wireframe(XX, YY, DensityMax2Dsorted) #, cmap=plt.cm.Blues) #FAILS: WHY do we have ground values ?! 
        #CS = ax.plot_surface(XX, YY, DensityMax2DsortedShape) #, cmap=plt.cm.Blues) #FAILS: WHY do we have ground values ?! 
        #CS = ax.plot_trisurf(BrutX, BrutY, BrutZ, cmap=plt.cm.Blues, antialiased=True, edgecolor='none') #WORKS
        #ax.view_init(90, 90)
        #See https://matplotlib.org/gallery/images_contours_and_fields/tricontour_smooth_user.html#sphx-glr-gallery-images-contours-and-fields-tricontour-smooth-user-py
        plt.triplot(triang, lw=0.1, color='black')
        levels = np.linspace(np.min(BrutZ), np.max(BrutZ), 10)
        #cmap = cm.get_cmaps(name='terrain', lut=None)
        cmap = cm.get_cmap(name='Blues', lut=None)
        plt.tricontourf(tri_refi, z_test_refi, levels=levels, cmap=cmap)
        #plt.tricontour(tri_refi, z_test_refi, levels=levels,
                #colors=['0.25', '0.5', '0.5', '0.5', '0.5'],
                #linewidths=[1.0, 0.5, 0.5, 0.5, 0.5])
        plt.colorbar()
        plt.tight_layout()
    except:
        #CS = -1
        print Header+"** Warning: failed to produce one of the plots."
    
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
    print Header+"** Info: saved "+filename+"eps|png."
    if(ShowFinalPlot==1):
      plt.show()
  
  del DensityMax2D, DensityMax2Dsorted
#}}}

##### 2D PLOTTING Density(wavelength1, wavelength2) for one single field value? 
#Mode 2: fixed fields, same pulse durations, then plot Ne_final(Color1, Color2)
def plotEffectOfDoublePulse_Color1_Color2_SingleFieldStrength(PlottingTable, FieldStrength, WhatToPlot, PulseDuration=20e-15, Tolerance=0.01, Egap_SI=2.56*e): #{{{
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
  
  #preparing the 2D plot
  DensityMax2D = np.empty((0,3))
  CurrentMax2D = np.empty((0,3))
  EnergyMax2D  = np.empty((0,3))
  Wpi2D        = np.empty((0,3))
  for Color1 in wavelength_uniq: #selecting color1
    PlottingDoublePulses_SelectedColor1 = FilterDatabase(PlottingDoublePulses, Color1, IndexTable['wavelengths1'])
    for Color2 in wavelength2_uniq: #selecting color2
      # Selecting Color 2
      PlottingDoublePulses_SelectedColor2 = FilterDatabase(PlottingDoublePulses_SelectedColor1, Color2, IndexTable['wavelengths2'])
      
      # Constructing the density table
      runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingDoublePulses_SelectedColor2)
      
      print "Color1: "+str(Color1)+" m, Color2: "+str(Color2)+" m."
      if(len(Ne_final) >= 1):
        ToBeAdded        = np.array([wavelengths1[0], wavelengths2[0], Ne_final[0]])
        ToBeAddedCurrent = np.array([wavelengths1[0], wavelengths2[0], TotalCurrent_final[0]])
        ToBeAddedEnergy  = np.array([wavelengths1[0], wavelengths2[0], Energy_max[0]])
        ToBeAddedWpi     = np.array([wavelengths1[0], wavelengths2[0], W_pi[0]])
        print ToBeAdded, ToBeAddedCurrent, ToBeAddedEnergy, ToBeAddedWpi
        DensityMax2D=np.vstack((DensityMax2D, ToBeAdded))
        CurrentMax2D=np.vstack((CurrentMax2D, ToBeAddedCurrent))
        EnergyMax2D =np.vstack((EnergyMax2D, ToBeAddedEnergy))
        Wpi2D       =np.vstack((Wpi2D, ToBeAddedWpi))
        del ToBeAdded, ToBeAddedCurrent, ToBeAddedEnergy, ToBeAddedWpi
        
  # Stack is built. Export, show and plot on a heat map. 
  print "** Info: exporting the sub-plot2D: Ne_final(wavelength, Field), TotalCurrent_Final(), Energy_max()..."
  print "** Info: useful to visualize available data with GNUplot. "
  filename         = "Bicolor-Field"+str(round(FieldStrength*1E-9))+"VperNm-NeFinal"
  filename_current = "Bicolor-Field"+str(round(FieldStrength*1E-9))+"VperNm-JeFinal"
  filename_energy  = "Bicolor-Field"+str(round(FieldStrength*1E-9))+"VperNm-Energy"
  filename_wpi     = "Bicolor-Field"+str(round(FieldStrength*1E-9))+"VperNm-Wpi"
  
  # Sorting the data before storing and plotting
  DensityMax2Dsorted = DensityMax2D[np.lexsort((DensityMax2D[:,2], DensityMax2D[:,1], DensityMax2D[:,0]))]
  CurrentMax2Dsorted = CurrentMax2D[np.lexsort((CurrentMax2D[:,2], CurrentMax2D[:,1], CurrentMax2D[:,0]))]
  EnergyMax2Dsorted  = EnergyMax2D[np.lexsort((EnergyMax2D[:,2], EnergyMax2D[:,1], EnergyMax2D[:,0]))]
  Wpi2Dsorted        = Wpi2D[np.lexsort((Wpi2D[:,2], Wpi2D[:,1], Wpi2D[:,0]))]
  #DensityMax2D.sort(axis=1) ##WARNING: THIS IS OUT OF CONTROL. DO NOT USE. 
  
  ExportToTxt(DensityMax2Dsorted, filename+".csv") #this is voluntarily unsorted to avoid damage
  ExportToTxt(CurrentMax2Dsorted, filename_current+".csv") #this is voluntarily unsorted to avoid damage
  ExportToTxt(EnergyMax2Dsorted,  filename_energy+".csv") #this is voluntarily unsorted to avoid damage
  ExportToTxt(Wpi2Dsorted,        filename_wpi+".csv") #this is voluntarily unsorted to avoid damage
  
  ## Lazy to build up a new function for that for now. 
  if(WhatToPlot == 1): 
      DensityMax2Dsorted = CurrentMax2Dsorted; 
      filename = filename_current
      QuantityTitle="Shift current $J_{end}(\lambda_1,\lambda_2)$"
  elif(WhatToPlot == 2): 
      DensityMax2Dsorted = EnergyMax2Dsorted
      filename = filename_energy
      QuantityTitle=r"Total energy $\xi_{end}(\lambda_1,\lambda_2)$"
  elif(WhatToPlot == 3): 
      DensityMax2Dsorted = Wpi2Dsorted
      filename = filename_wpi
      QuantityTitle=r"Excitation rate $w_{PI}(\lambda_1,\lambda_2)$"
  else: 
      QuantityTitle=r"Excited density $N_{exc}(\lambda_1,\lambda_2)$"
  
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
    plt.title(QuantityTitle+', '+r'$E_{1,2}=$'+str(FieldStrength*1E-9)+' V/nm')
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

## Prepare all plots for single-color excitation with VP Zhukov [Phys. Rev. B]
def plot_Zhukov_SingleColor_LDA(): 
    Egap_SI = 2.56*e
    PulseType  = 3 #STH: 2, StH: 3
    wavelength = 800e-9
    ShowPlot   = 0
    PlotWpiAxis = True

    ## Plots effect of the STH Rise Time. We even found an incomplete computations thanks to this. Works super well, then. 
    RiseCycles = 0 #0: undefined, mixes everything. Basically, we must redefine all runs which had no well defined RiseCycle so far. 
    Efield     = 1E9
    plotEffectOfRiseTime(PlottingTable, wavelength, 20e-15, Efield, PulseType, ShowPlot, Egap_SI)
    
    ## Effects of CEP on single pulse excitation, current and energy, and harmonic spectrum
    RiseCycles = 2 #NOTE: a number different from 0 must be chosen before publication. 
    Efield     = 1E9
    plotEffectOfCEP(PlottingTable, wavelength, 20e-15, Efield, PulseType, RiseCycles, ShowPlot, Egap_SI)
    #plotEffectOfCEP(PlottingTable, wavelength, 10e-15, Efield, PulseType, RiseCycles, ShowPlot)
    Efield     = 5E9 #NOTE: Saturation may prevent to see any difference. 
    plotEffectOfCEP(PlottingTable, wavelength, 20e-15, Efield, PulseType, RiseCycles, ShowPlot, Egap_SI)
    #plotEffectOfCEP(PlottingTable, wavelength, 10e-15, Efield, PulseType, RiseCycles, ShowPlot)
    
    #NOTE: === LOCAL FILTERING
    #PlottingTable_filtered = FilterDatabaseRemove(PlottingTable, 30e-15, IndexTable['pulses1'])
    PlottingTable_filtered = PlottingTable
    
    ## Effect of intensity figures
    plotEffectOfIntensity_CEP(PlottingTable_filtered, wavelength, 20e-15, False, 0.0, PulseType, RiseCycles, ShowPlot, Egap_SI)
    #plotEffectOfIntensity_CEP(PlottingTable_filtered, wavelength, 10e-15, False, 0.0, PulseType, RiseCycles, ShowPlot)
    #plotEffectOfIntensity_CEP(PlottingTable_filtered, wavelength, 30e-15, False, 0.0, PulseType, RiseCycles, ShowPlot)
    
    plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI) #This one construct the ionization probabilities! 
    wavelength = 483e-9 #Resonant wavelength results for LDA.
    #wavelength = 407e-9 #Resonant wavelength results for TB09. 
    plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI) #This one construct the ionization probabilities! 
    wavelength = 484e-9 #Resonant wavelength results for LDA.
    plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI) #This one construct the ionization probabilities! 
    
    ## READY
    wavelength = 800e-9
    plotEffectOfIntensity_wavelengths(PlottingTable_filtered, 20e-15, 0., False, wavelength, PulseType, RiseCycles, ShowPlot, Egap_SI) 
    #plotEffectOfIntensity_wavelengths(PlottingTable_filtered, 10e-15, 0., False, wavelength, PulseType, RiseCycles, ShowPlot)
    plotEffectOfIntensity_wavelengths(PlottingTable_filtered, 30e-15, 0., False, wavelength, PulseType, RiseCycles, ShowPlot)
    
    ShowPlot=1
    #plotEffectOfWavelength_2dMap(PlottingTable, 20e-15, PulseType, RiseCycles, ShowPlot)
    plotEffectOfWavelength_Tolerance(PlottingTable, 0.03, 20e-15, 0., PulseType, RiseCycles, ShowPlot, Egap_SI)
    ### OLD STUFFS
    #plotEffectOfIntensity_CEP(PlottingTable, wavelength, 10e-15, PulseType, RiseCycles, ShowPlot)
    ## Plots the N_exc(wavelength, field) (2D plot) #NOTE: this can be compared with GrojoPRB2013. 

def plot_Zhukov_SingleColor_TB09():
    Egap_SI = 3.04*e
    PulseType  = 3 #STH: 2, StH: 3
    wavelength = 800e-9
    ShowPlot   = 0
    PlotWpiAxis = True
    ## Plots effect of the STH Rise Time. We even found an incomplete computations thanks to this. Works super well, then. 
    #RiseCycles = 0 #0: undefined, mixes everything. Basically, we must redefine all runs which had no well defined RiseCycle so far. 
    #Efield     = 1E9
    #plotEffectOfRiseTime(PlottingTable, wavelength, 20e-15, Efield, PulseType, ShowPlot)
    
    ## Effects of CEP on single pulse excitation, current and energy, and harmonic spectrum
    RiseCycles = 2 #NOTE: a number different from 0 must be chosen before publication. 
    #Efield     = 1E9
    #plotEffectOfCEP(PlottingTable, wavelength, 20e-15, Efield, PulseType, RiseCycles, ShowPlot)
    ##plotEffectOfCEP(PlottingTable, wavelength, 10e-15, Efield, PulseType, RiseCycles, ShowPlot)
    #Efield     = 5E9 #NOTE: Saturation may prevent to see any difference. 
    #plotEffectOfCEP(PlottingTable, wavelength, 20e-15, Efield, PulseType, RiseCycles, ShowPlot)
    ##plotEffectOfCEP(PlottingTable, wavelength, 10e-15, Efield, PulseType, RiseCycles, ShowPlot)
    
    #NOTE: === LOCAL FILTERING
    #PlottingTable_filtered = FilterDatabaseRemove(PlottingTable, 30e-15, IndexTable['pulses1'])
    PlottingTable_filtered = PlottingTable
    
    ## Effect of intensity figures
    #plotEffectOfIntensity_CEP(PlottingTable_filtered, wavelength, 20e-15, False, 0.0, PulseType, RiseCycles, ShowPlot)
    #plotEffectOfIntensity_CEP(PlottingTable_filtered, wavelength, 10e-15, False, 0.0, PulseType, RiseCycles, ShowPlot)
    #plotEffectOfIntensity_CEP(PlottingTable_filtered, wavelength, 30e-15, False, 0.0, PulseType, RiseCycles, ShowPlot)
    
    # Compute ionization probability, effect of the pulse duration
    ShowPlot=1
    plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI)
    
    #This one construct the ionization probabilities! 
    #wavelength = 483e-9 #Resonant wavelength results for LDA.
    wavelength = 407.58e-9 #Resonant wavelength results for TB09. 
    plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI) #This one construct the ionization probabilities! 
    exit()
    #wavelength = 484e-9 #Resonant wavelength results for LDA.
    #plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot) #This one construct the ionization probabilities! 
    ## READY
    wavelength = 800e-9
    plotEffectOfIntensity_wavelengths(PlottingTable_filtered, 20e-15, 0., PlotWpiAxis, wavelength, PulseType, RiseCycles, ShowPlot, Egap_SI) 
    #plotEffectOfIntensity_wavelengths(PlottingTable_filtered, 10e-15, 0., PlotWpiAxis, wavelength, PulseType, RiseCycles, ShowPlot)
    #plotEffectOfIntensity_wavelengths(PlottingTable_filtered, 30e-15, 0., PlotWpiAxis, wavelength, PulseType, RiseCycles, ShowPlot)
    
    ShowPlot=1
    #plotEffectOfWavelength_2dMap(PlottingTable, 20e-15, PulseType, RiseCycles, ShowPlot)
    plotEffectOfWavelength_Tolerance(PlottingTable, 0.03, 20e-15, 0., PulseType, RiseCycles, ShowPlot, Egap_SI)
    ### OLD STUFFS
    #plotEffectOfIntensity_CEP(PlottingTable, wavelength, 10e-15, PulseType, RiseCycles, ShowPlot)
    ## Plots the N_exc(wavelength, field) (2D plot) #NOTE: this can be compared with GrojoPRB2013. 



## Prepare all plots for bicolor paper with V. Zhukov [1 Phys. Rev. Lett + 1 Phys. Rev. A]
def plot_Zhukov_BicolorMaps():
    for tau in [20e-15, 30e-15]:
        PulseType  = 3; RiseCycles = 2 #StH corrected. 
        #PulseType  = 2;  RiseCycles = 0 #bugged STH. 
        
        ## Commensurable
        plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, 800e-9, 1600e-9, tau, PulseType, RiseCycles,  0)
        plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, 800e-9, 1600e-9, tau, PulseType, RiseCycles,  1)
        plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, 800e-9, 1600e-9, tau, PulseType, RiseCycles,  2)
        plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, 800e-9, 1600e-9, tau, PulseType, RiseCycles,  3)
        plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, 800e-9, 1600e-9, tau, PulseType, RiseCycles,  4)
        ## Incommensurable                                                                                  
        plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, 800e-9, 1030e-9, tau, PulseType, RiseCycles,  0)
        plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, 800e-9, 1030e-9, tau, PulseType, RiseCycles,  1)
        plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, 800e-9, 1030e-9, tau, PulseType, RiseCycles,  2)
        plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, 800e-9, 1030e-9, tau, PulseType, RiseCycles,  3)
        plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, 800e-9, 1030e-9, tau, PulseType, RiseCycles,  4)

## Prepare all plots for shift current paper with A. Rubio
def plot_ShiftCurrents():
    for tau in [30e-15, 50e-15]:
        PulseType  = 1  #SinSquared
        RiseCycles = 0  #We dont care  
        wavelength1 = 3100e-9
        for wavelength2 in [1.55e-6]:
            ## Commensurable
            plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  0)
            plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  1)
            plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  2)
            plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  3)
            ## Incommensurable                                                                                        
            plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  0)
            plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  1)
            plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  2)
            plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  3)

    
    
    

# VALID but some degeneracies !
#plotEffectOfIntensity(PlottingTable, 10e-15, 1, 1, 1) #validated
#plotEffectOfIntensity_wavelengths(PlottingTable, 20e-15, 1, 0) #validated, degeneracies are due to undefined RiseTime
#plotEffectOfIntensity_wavelengths(PlottingTable, 30e-15, 0) #validated
#plotEffectOfIntensity_wavelengths(PlottingTable, 40e-15, 0) #validated
#plotEffectOfIntensity_wavelengths(PlottingTable, 50e-15, 0) #validated
#plotEffectOfIntensity_PulseDuration(PlottingTable, 800e-9, 1) #TODO: remove the different CEPs and keep only tau dependence. 


# TYPICAL RUNS
#plotEffectOfWavelength_2dMap(PlottingTable) #VALID works, but not enough data for now. 
#plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple
#plotEffectOfDoublePulse_Color1_Color2_SingleFieldStrength(PlottingTable, 5E9, 0) #VALID, but don't apply to TOO LARGE wavelength sets
#plotEffectOfDoublePulse_Color1_Color2_SingleFieldStrength(PlottingTable, 5E9, 1) #VALID, but don't apply to TOO LARGE wavelength sets
#plotEffectOfDoublePulse_Color1_Color2_SingleFieldStrength(PlottingTable, 5E9, 2) #VALID, but don't apply to TOO LARGE wavelength sets
#plotEffectOfDoublePulse_Color1_Color2_SingleFieldStrength(PlottingTable, 5E9, 3) #VALID, but don't apply to TOO LARGE wavelength sets

# DEVELOPMENT RUNS
#plotEffectOfWavelength_Tolerance(PlottingTable) #considers a tolerance to plot various fields
#plotEffectOfCEP(PlottingTable)
# plotEffectOfFrozenXC(PlottingTable)

#plotEffectOfDoublePulse_Efield1_Efield2_batchAllColors(PlottingTable) #BUGGY
#plotEffectOfDoublePulse_sameAmplitude_samePulseDurations(PlottingTable)
#plotEffectOfDoublePulse_diffWavelengths(PlottingTable)


#====================================== EXECUTION ZONE ============================
## Zhukov single color plotting [Phys. Rev. B]
#plot_Zhukov_SingleColor_LDA()
plot_Zhukov_SingleColor_TB09()

## Zhukov bicolor plotting
#plot_Zhukov_BicolorMaps()

## Rubio shift current plotting
#plot_ShiftCurrents()

print Header+"End of program."
exit()

