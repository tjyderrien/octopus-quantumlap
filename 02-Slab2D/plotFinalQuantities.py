#!/usr/bin/env python2
#-*- coding: utf-8 -*-

# Copyright (C) 2016-2020 T. J.-Y. Derrien
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

## @package plotFinalQuantities
# Collects the generated time-dependent octopus simulations results from all subfolders (can last minutes), builds a table as function of the laser parameters, final densities, final current density, and final energies after simulation. Then prepares a bunch of scientific results for publications with one and two color simulations. 

## IMPORT LIBRARIES
import numpy as np
from numpy import genfromtxt, loadtxt, chararray
#from scipy import interpolate
import matplotlib.pyplot as pplt
import matplotlib.pylab  as pl
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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#from matplotlib.tri import Triangulation, UniformTriRefiner
from scipy.interpolate import interp2d
from scipy.optimize import curve_fit
from colorama import Fore
from colorama import Style

from math import ceil

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
rc('text', usetex=False) #True: Does not work on Draco. Just put False, then, and it will be smooth. 
mp.rcParams['legend.numpoints'] = 1

Header="[plotFinalQuantities] "

EnablePeakFinder           = 0
ShowFinalPlot              = 0

BicolorE1E2plotDictionnary={'Nexc': 0, 'E2': 4, 'J': 1, 'Energy': 2, 'Wpi': 3, 'NexcGrad': 5}

Density=5E28 #specific to Si #TODO: compute automatically? 
n0 = 4.*Density #maximum density of electrons in Si
xscale = 1E-4 #convert W/m2 to W/cm2. 

#Egap_SI = 2.56*e #LDA: 2.56 eV, TB09: 3.04 eV. [Waroquiers and Gonze, Phys. Rev. B 87, 075121 (2013)] #TODO: this will depend on which functional we have used! TB09 or LDA ! 

# Capture the data
# For histocal reasons, two locations are possilbe.
path1='*/*/output_iter/' #valid for the runs computed from a same k-grid folder. 
path2='*/*/td.general/'

if(EnablePeakFinder): 
  from scipy.signal import find_peaks

PlottingTable = libExcitationDatabase.GetExcitationDatabase(path2) #completes the existing database, or fetches the new one in all the subfolders where this is executed. 

print PlottingTable

PlottingTable = FilterDatabase(PlottingTable, 0, IndexTable['FieldInduced']) #removes all induced fields data

#: TODO one has to extract the volume! But how to proceed when it is time-dependent? 

runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingTable)

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
SizeX = 12.
#SizeY = SizeX / GoldNumber
SizeY = 6.
OneSizeFitsAll="small"  #caption size for all legends of physics or rocknroll

## Plot Nexc(RiseTime; single tau, 1 curve per Wavelength)
def plotEffectOfRiseTime(PlottingTable, pulseWavelength=800e-9, pulseDuration=20e-15, Efield1=5E9, PulseShape=3, ShowFinalPlot=0, Egap_SI=2.56*e): #{{{
  print "#====== Plot: SINGLE PULSE, EFFECT OF THE RISE TIME, 1 WAVELENGTH, 1 POLARIZATIONS: final density as function of RISETIME ======"
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
  pplt.title(r'$n_{\mathrm{exc}}(RiseCycles)$, $\tau=$'+str(int(pulseDuration*1E15))+' fs.')
  # Number of sweeped parameters: 3 (Efield, Delay, CrossPolarization)
  #pplt.scatter(FieldToIntensity(Efield1), Ne_final, c=Delay, s=(30.+CrossPolarization*50.), label=r'$N_{exc}$[Field]')
  plot1 = [] #np.empty((0,0))
  for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
    PlottingSinglePulses_CEP = FilterDatabase(PlottingSinglePulses, CEP1_value, IndexTable['CEP1']) #select the local CEP
    #TODO: filter by field strength ! Maybe take a piece of code from tolerance-based selector in another function. #TODO: Use subfunctions otherwise we will get crazy. 
    ## prints the current state of selection and stores the values of sweeped parameters
    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_CEP)
    
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
  plt.ylabel(r'Final density (m$^{-3}$)')
  #  pplt.ylim((1E16,4.*5e28))
  #ax1.set_xscale('log')
  #ax1.set_yscale('log')
  #plot1_list = plot1.tolist()
    #print Efield1, CEP1_value
  
  ax12 = ax1.twinx()
  for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
    PlottingSinglePulses_CEP = FilterDatabase(PlottingSinglePulses, CEP1_value, IndexTable['CEP1']) #select the local CEP
    
    ## prints the current state of selection and stores the values of sweeped parameters
    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_CEP)
    
    #plot_field, = ax1.plot(FieldToIntensity(Efield1), Ne_final, '-', label=r'$N_{exc}^{final}$, $\phi_{CE}$='+str(round(CEP1_value*180./pi)))
    #lists = sorted(itertools.izip(*[FieldToIntensity(Efield1), Ne_final]))
    lists = sorted(itertools.izip(*[RiseTime, Ne_max, Ne_final, pulses1])) #NOTE: Ne_max makes more sense than Ne_final
    RiseTime_sorted, Ne_max_sorted, Ne_final_sorted, pulses1_sorted = list(itertools.izip(*lists))
    plt.plot(RiseTime_sorted, np.divide(Ne_final_sorted, pulses1_sorted), 'o') #, label=r'$N_{exc}/\tau$') #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
  
  #ax12.set_yscale('log')
  ax12.set_ylabel(r'$w_{\mathrm{PI}}$ (m$^{-3}$.s$^{-1}$)')

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

## Plots Nexc(intensity; single tau, 1 curve per wavelength) | Wpi(intensity)
def plotEffectOfIntensity_wavelengths(PlottingTable, pulseDuration=20e-15, CEP=0., RestrictWavelength1=False, wavelength_value=800e-9, PulseType=3, RiseTime=2, ShowFinalPlot=0, Egap_SI = 2.56*e): #{{{
  PlotVirtualExcitation = False
  print "#====== Plot: SINGLE PULSE, VARIOUS WAVELENGTHS, SAME PULSE DURATION, SAME POLARIZATIONS: final density as function of field ====== "
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
  #if(len(wavelength_uniq) == 1): #single color simulations
  # Complexity level: 1
  plot1=[]
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

    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_wavelength)
    lists = sorted(itertools.izip(*[wavelengths1,FieldToIntensity(Efield1), Ne_final, pulses1]))
    wavelengths_sorted, Intensity_sorted, Ne_final_sorted, pulses1_sorted = list(itertools.izip(*lists))
    plotNe, = pplt.plot(Intensity_sorted, Ne_final_sorted, 'o', label=r'$N_{exc}$, $\lambda=$'+str(round(wavelength_value*1E9,0))+' nm')
  #}}}
  #ax1.legend(loc='best', fancybox=True, framealpha=0) #bbox_to_anchor=(1.05, 1), 
  if(not PlotVirtualExcitation): 
    pplt.xlabel(r'Peak intensity (W/cm$^2$)')
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

    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_wavelength)

    lists = sorted(itertools.izip(*[wavelengths1,FieldToIntensity(Efield1), Ne_final, pulses1]))
    wavelengths_sorted, Intensity_sorted, Ne_final_sorted, pulses1_sorted = list(itertools.izip(*lists))

    pplt.plot(Intensity_sorted, np.divide(Ne_final_sorted, (np.multiply(2.,pulses1_sorted))), 'o') #, label=r'$N_{exc}/\tau$, $\lambda=$'+str(round(wavelength_value*1E9))+' nm')
  #}}}
  ax12.set_yscale('log')
  pplt.ylabel('$w_{\mathrm{PI}}$ (m$^{-3}$.s$^{-1}$)')
  
  # Second panel of the figure
  if(PlotVirtualExcitation):
    ax2 = pplt.subplot(plotB)
    #pplt.title('rVirtual exc. $\zeta$')
    pplt.plot(FieldToIntensity(Efield1), irreversibility, 'ob', label='Virtual exc.')
    pplt.xlabel(r'Peak intensity (W/cm$^2$)')
    pplt.ylabel(r'Virtual exc. $\zeta$')
    pplt.ylim((0,1))
    ax2.set_xscale('log')
    ax22 = ax2.twinx()
    for wavelength_value in wavelength_uniq:
      # select the corresponding wavelength in results
      PlottingSinglePulses_wavelength = FilterDatabase(PlottingSinglePulses, wavelength_value, IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts.

      runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_wavelength)
      pplt.plot(FieldToIntensity(Efield1), gammaKeldysh(Egap_SI, meff, Efield1, wavelengths1), 'or', label=r'Adiabaticity, $\lambda=$'+str(round(wavelength_value*1E9, 0))+' nm')
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

## Plots (N_exc | Wpi | Energy) (peak intensity), (1 curve per pulse duration)
# Nexc(I) is used to extract 1 to 7 photon absorption rates, given in caption. 
# Wpi is compared to Keldysh analytical solid theory, and numerically integrated Keldysh if available. # Wpi is also compared to tunneling limit case. 
# ElectronEnergy (intensity) is given. 
def plotEffectOfIntensity_PulseDuration(PlottingTable, pulseWavelength=800e-9, CEP0=0., PulseType=3, RiseCycles=2, PlotWpiAxis=False, ShowFinalPlot=0, Egap_SI=2.56*e, Egap_SI_indirect=0.51*e, OpticalIndex=1.0, MoveIon=0): #{{{
  PlotVirtualExcitation = False #Not advised! This is not a real observable. 
  PlotEnergy = True #PRL: False. #PRB: True.
  EnableKeldysh = True
  EnableStarkEffect = False #Stark effect on tunneling? Disabled in original Keldysh, Gruzdev and Kaiser theories. 
  ShowPulseDurations=False
  
  Ntotal=5.0E28 
  meff_Waroquiers_LDA = 0.2226 #0.2226 #Waroquiers, D. et al. Physical Review B, 2013, 87, 075121.
  meff = meff_Waroquiers_LDA #0.14 #used to adjust the w_PI model 
  meff_Tunnel = meff_Waroquiers_LDA #fits well if using INDIRECT band gap. 
  
  
  print "====== Plot: SINGLE PULSE, SAME WAVELENGTHS, SAME POLARIZATIONS: final density as function of field and pulse duration ======"
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
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, MoveIon, IndexTable['MoveIon']) #filters the enabled or disabled movement of ions. 
  
  if(pulseWavelength==1600e-9): #NOTE: enable if you wish to remove 20 fs results for 1600 nm.
      PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, 21e-15, IndexTable['pulses1'])
  if(pulseWavelength==800e-9): #NOTE: enable if you wish to remove 10 fs results for 800 nm.
      PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, 11e-15, IndexTable['pulses1'])
  
  ## NOTE: For clarity, I will filter results that are supplementary for the paper's figures. 
  ## Let's remove the results of 1600 nm with 20 fs pulse duration. We can keep 30 fs. 

  ## prints the current state of selection
  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = CountTheNumberOfLaserParameters(PlottingSinglePulses)

  # Simple pulse, single color
  #plotA = 311; plotB = 312; plotC = 313 #Horizontal
  if(not PlotEnergy):
      plotA=121; plotB=122;
      pplt.figure(figsize=(SizeX/3.*2.,SizeY))
  else:
      plotA = 131; plotB = 132; plotC = 133 #Vertical
      pplt.figure(figsize=(SizeX,SizeY))
      
  #plotA = 221; plotB = 224; plotC = 223; plotD = 222 # 2x2 for huge caption of ax1
  
  EnergyMin = 0e0 #-1E-7
  NumberOfPhotons_round = 99 #we need to initialize it here
  
  ax1 = pplt.subplot(plotA)
  
  plot1 = []
  Symbols={10e-15: '+', 20e-15: 'o', 30e-15: '^'}
  Colors={10e-15: 'g', 20e-15: 'r', 30e-15: 'b'}
  for pulse_value in pulses_uniq: #{{{
    # select the corresponding wavelength in results
    PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses, pulse_value, IndexTable['pulses1'])
    PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses_pulses, pulse_value, IndexTable['pulses2']) #keeps only pulses of the same duration.
    
    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIonRun, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_pulses)
    
    print Header+"** Info: saving the data in CSV file for fitting. "
    filename = "EffectOfPulseIntensity-Duration-"+str(pulseWavelength*1E9)+"nm-tau"+str(int(pulse_value*1E15))+"fs"
    print Header+"** Info: saving in "+filename+".csv"
    ExportToTxt(PlottingSinglePulses_pulses, filename+".csv")
    lists = sorted(itertools.izip(*[pulses1, Efield1, FieldToIntensity(Efield1), Ne_final]))
    pulses1_sorted, Efield1_sorted, Intensity_sorted, Ne_final_sorted = list(itertools.izip(*lists))
    #print Header+"Ready to plot: "
    
    plot_field, = ax1.plot(np.multiply(xscale,Intensity_sorted), Ne_final_sorted, Colors[pulse_value]+Symbols[pulse_value], label=r'$n_{\mathrm{exc}}^{\mathrm{TDDFT}}$, $\tau=$'+str(int(pulse_value*1E15))+' fs')
    plot1 += [plot_field]
    
    #==== CURVE FITTING SHOULD COME HERE === as it is for each pulse duration. 
    print Header+"** Info: Computation of multi-photon excitation"
    
    ### Multiphotonic absorption law (used for fitting td-dft results)
    # d e f MultiPhotonAbsorption(intensity, sigma_n):
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
    
    ## 4-photon absorption fitting function
    def FourPhotonAbsorption(intensity, sigma_n): 
    #{{{
      n=4; #wavelength=800e-9; tau_FWHM=20e-15; 
      omega = 2.*np.pi*c/wavelengths1[0]
      tau_FWHM = pulses1_sorted[0]
      #nexc_mpi = MultiPhotonAbsorption(intensity, sigma_n, n, wavelength, tau_FWHM)
      nexc_mpi = sigma_n * intensity**n / (n*hbar*omega) * tau_FWHM
      return nexc_mpi
    #}}}
    
    ## 5-photon absorption fitting function
    def FivePhotonAbsorption(intensity, sigma_n): 
    #{{{
      n=5; #wavelength=800e-9; tau_FWHM=20e-15; 
      omega = 2.*np.pi*c/wavelengths1[0]
      tau_FWHM = pulses1_sorted[0]
      #nexc_mpi = MultiPhotonAbsorption(intensity, sigma_n, n, wavelength, tau_FWHM)
      nexc_mpi = sigma_n * intensity**n / (n*hbar*omega) * tau_FWHM
      return nexc_mpi
    #}}}
    
    ## 6-photon absorption fitting function
    def SixPhotonAbsorption(intensity, sigma_n): 
    #{{{
      n=6; #wavelength=800e-9; tau_FWHM=20e-15; 
      omega = 2.*np.pi*c/wavelengths1[0]
      tau_FWHM = pulses1_sorted[0]
      #nexc_mpi = MultiPhotonAbsorption(intensity, sigma_n, n, wavelength, tau_FWHM)
      nexc_mpi = sigma_n * intensity**n / (n*hbar*omega) * tau_FWHM
      return nexc_mpi
    #}}}
    
    ## 7-photon absorption fitting function
    def SixPhotonAbsorption(intensity, sigma_n): 
    #{{{
      n=7; #wavelength=800e-9; tau_FWHM=20e-15; 
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
    def FitTwoPhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, n_photons, TwoPhotonIntensityRegime_min, TwoPhotonIntensityRegime_max, plotMe=True): #TODO: pass function name in argument
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
        if(ShowPulseDurations):
            tau_text=r', $\tau=$'+str(int(pulse_value*1E15))+' fs'
        else: 
            tau_text=""
        if(plotMe):
            plot_field, = ax1.plot(np.multiply(xscale,Intensity_sorted_mask), Nexc_2PI, Colors[pulse_value]+'--', label=r'$\sigma_2=$'+str('%.4g'%TwoPhotonAbsorption_params[0][0])+tau_text) #NOTE add value of sigma_2 here
        else: 
            plot_field = []
        del Nexc_2PI, TwoPhotonAbsorption_params, Intensity_sorted_mask
        return plot_field
    #}}}
    
    ## Filters a part of the data and fit using the given function
    def FitOnePhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, n_photons, OnePhotonIntensityRegime_min, OnePhotonIntensityRegime_max, plotMe=True): #TODO: pass function name in argument
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
        if(ShowPulseDurations):
            tau_text=r', $\tau=$'+str(int(pulse_value*1E15))+' fs'
        else: 
            tau_text=""
        if(plotMe):
            plot_field, = ax1.plot(np.multiply(xscale,Intensity_sorted_mask), Nexc_1PI, Colors[pulse_value]+'-.', label=r'$\sigma_1=$'+str('%.4g'%OnePhotonAbsorption_params[0][0])+tau_text) #NOTE add value of sigma_1 here
        else:
            plot_field = []
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
        
        if(ShowPulseDurations):
            tau_text=r', $\tau=$'+str(int(pulse_value*1E15))+' fs'
        else: 
            tau_text=""
            
        print Header+"** Info: adding the fitted curve to the plot"
        Nexc_3PI = ThreePhotonAbsorption(Intensity_sorted_mask, ThreePhotonAbsorption_params[0])
        #print Nexc_2PI
        plot_field, = ax1.plot(np.multiply(xscale,Intensity_sorted_mask), Nexc_3PI, Colors[pulse_value]+':', label=r'$\sigma_3=$'+str('%.4g'%ThreePhotonAbsorption_params[0][0])+tau_text) #NOTE adds value of sigma_3 here
        del Nexc_3PI, ThreePhotonAbsorption_params, Intensity_sorted_mask
        return plot_field
    #}}}
    
    def FitFourPhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, n_photons, FourPhotonIntensityRegime_min, FourPhotonIntensityRegime_max): #TODO: pass function name in argument instead of copying this function... :-(
    #{{{
        print Header+"** Info: Filtering the fitting region..."
        Intensity_sorted_np = np.array(Intensity_sorted) #filtering works for NumPy array, not for lists. 
        #condition = (Intensity_sorted_np > TwoPhotonIntensityRegime_min) and (Intensity_sorted_np < TwoPhotonIntensityRegime_max)
        Intensity_sorted_mask = Intensity_sorted_np[(Intensity_sorted_np > FourPhotonIntensityRegime_min) & (Intensity_sorted_np < FourPhotonIntensityRegime_max)]
        
        Ne_final_sorted_np  = np.array(Ne_final_sorted)
        Ne_final_sorted_mask  = Ne_final_sorted_np[(Intensity_sorted_np > FourPhotonIntensityRegime_min) & (Intensity_sorted_np < FourPhotonIntensityRegime_max)]
        print Intensity_sorted_mask, Ne_final_sorted_mask #Works well
        
        FourPhotonAbsorption_params = curve_fit(FourPhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
        #ThreePhotonAbsorption_params = curve_fit(MultiPhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
        print Header+"** Ionization parameters for tau="+str(int(1e15*pulse_value))+" fs: "+str(FourPhotonAbsorption_params[0])
        
        if(ShowPulseDurations):
            tau_text=r', $\tau=$'+str(int(pulse_value*1E15))+' fs'
        else: 
            tau_text=""
            
        print Header+"** Info: adding the fitted curve to the plot"
        Nexc_4PI = FourPhotonAbsorption(Intensity_sorted_mask, FourPhotonAbsorption_params[0])
        #print Nexc_2PI
        plot_field, = ax1.plot(np.multiply(xscale,Intensity_sorted_mask), Nexc_4PI, Colors[pulse_value]+'--', label=r'$\sigma_4=$'+str('%.4g'%FourPhotonAbsorption_params[0][0])+tau_text) #NOTE adds value of sigma_3 here
        del Nexc_4PI, FourPhotonAbsorption_params, Intensity_sorted_mask
        return plot_field
    #}}}
    
    def FitFivePhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, n_photons, FivePhotonIntensityRegime_min, FivePhotonIntensityRegime_max): #TODO: pass function name in argument instead of copying this function... :-(
    #{{{
        print Header+"** Info: Filtering the fitting region..."
        Intensity_sorted_np = np.array(Intensity_sorted) #filtering works for NumPy array, not for lists. 
        #condition = (Intensity_sorted_np > TwoPhotonIntensityRegime_min) and (Intensity_sorted_np < TwoPhotonIntensityRegime_max)
        Intensity_sorted_mask = Intensity_sorted_np[(Intensity_sorted_np > FivePhotonIntensityRegime_min) & (Intensity_sorted_np < FivePhotonIntensityRegime_max)]
        
        Ne_final_sorted_np  = np.array(Ne_final_sorted)
        Ne_final_sorted_mask  = Ne_final_sorted_np[(Intensity_sorted_np > FivePhotonIntensityRegime_min) & (Intensity_sorted_np < FivePhotonIntensityRegime_max)]
        print Intensity_sorted_mask, Ne_final_sorted_mask #Works well
        
        FivePhotonAbsorption_params = curve_fit(FivePhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
        #ThreePhotonAbsorption_params = curve_fit(MultiPhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
        print Header+"** Ionization parameters for tau="+str(int(1e15*pulse_value))+" fs: "+str(FivePhotonAbsorption_params[0])
        
        print Header+"** Info: adding the fitted curve to the plot"
        Nexc_5PI = FivePhotonAbsorption(Intensity_sorted_mask, FivePhotonAbsorption_params[0])
        #print Nexc_2PI
        
        if(ShowPulseDurations):
            tau_text=r', $\tau=$'+str(int(pulse_value*1E15))+' fs'
        else: 
            tau_text=""
        
        plot_field, = ax1.plot(np.multiply(xscale,Intensity_sorted_mask), Nexc_5PI, Colors[pulse_value]+'-.', label=r'$\sigma_5=$'+str('%.4g'%FivePhotonAbsorption_params[0][0])+tau_text) #NOTE adds value of sigma_3 here
        del Nexc_5PI, FivePhotonAbsorption_params, Intensity_sorted_mask
        return plot_field
    #}}}
    
    def FitSixPhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, n_photons, SixPhotonIntensityRegime_min, SixPhotonIntensityRegime_max): #TODO: pass function name in argument instead of copying this function... :-(
    #{{{
        print Header+"** Info: Filtering the fitting region..."
        Intensity_sorted_np = np.array(Intensity_sorted) #filtering works for NumPy array, not for lists. 
        #condition = (Intensity_sorted_np > TwoPhotonIntensityRegime_min) and (Intensity_sorted_np < TwoPhotonIntensityRegime_max)
        Intensity_sorted_mask = Intensity_sorted_np[(Intensity_sorted_np > SixPhotonIntensityRegime_min) & (Intensity_sorted_np < SixPhotonIntensityRegime_max)]
        
        Ne_final_sorted_np  = np.array(Ne_final_sorted)
        Ne_final_sorted_mask  = Ne_final_sorted_np[(Intensity_sorted_np > SixPhotonIntensityRegime_min) & (Intensity_sorted_np < SixPhotonIntensityRegime_max)]
        print Intensity_sorted_mask, Ne_final_sorted_mask #Works well
        
        SixPhotonAbsorption_params = curve_fit(SixPhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
        #ThreePhotonAbsorption_params = curve_fit(MultiPhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
        print Header+"** Ionization parameters for tau="+str(int(1e15*pulse_value))+" fs: "+str(SixPhotonAbsorption_params[0])
        
        print Header+"** Info: adding the fitted curve to the plot"
        Nexc_6PI = SixPhotonAbsorption(Intensity_sorted_mask, SixPhotonAbsorption_params[0])
        #print Nexc_2PI
        
        if(ShowPulseDurations):
            tau_text=r', $\tau=$'+str(int(pulse_value*1E15))+' fs'
        else: 
            tau_text=""
            
        plot_field, = ax1.plot(np.multiply(xscale,Intensity_sorted_mask), Nexc_6PI, Colors[pulse_value]+':', label=r'$\sigma_6=$'+str('%.4g'%SixPhotonAbsorption_params[0][0])+tau_text) #NOTE adds value of sigma_3 here
        del Nexc_6PI, SixPhotonAbsorption_params, Intensity_sorted_mask
        return plot_field
    #}}}
    
    def FitSevenPhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, n_photons, SevenPhotonIntensityRegime_min, SevenPhotonIntensityRegime_max): #TODO: pass function name in argument instead of copying this function... :-(
    #{{{
        print Header+"** Info: Filtering the fitting region..."
        Intensity_sorted_np = np.array(Intensity_sorted) #filtering works for NumPy array, not for lists. 
        #condition = (Intensity_sorted_np > TwoPhotonIntensityRegime_min) and (Intensity_sorted_np < TwoPhotonIntensityRegime_max)
        Intensity_sorted_mask = Intensity_sorted_np[(Intensity_sorted_np > SevenPhotonIntensityRegime_min) & (Intensity_sorted_np < SevenPhotonIntensityRegime_max)]
        
        Ne_final_sorted_np  = np.array(Ne_final_sorted)
        Ne_final_sorted_mask  = Ne_final_sorted_np[(Intensity_sorted_np > SevenPhotonIntensityRegime_min) & (Intensity_sorted_np < SevenPhotonIntensityRegime_max)]
        print Intensity_sorted_mask, Ne_final_sorted_mask #Works well
        
        SevenPhotonAbsorption_params = curve_fit(SevenPhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
        #ThreePhotonAbsorption_params = curve_fit(MultiPhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
        print Header+"** Ionization parameters for tau="+str(int(1e15*pulse_value))+" fs: "+str(SevenPhotonAbsorption_params[0])
        
        print Header+"** Info: adding the fitted curve to the plot"
        Nexc_7PI = SevenPhotonAbsorption(Intensity_sorted_mask, SevenPhotonAbsorption_params[0])
        #print Nexc_2PI
        
        if(ShowPulseDurations):
            tau_text=r', $\tau=$'+str(int(pulse_value*1E15))+' fs'
        else: 
            tau_text=""
            
        plot_field, = ax1.plot(np.multiply(xscale,Intensity_sorted_mask), Nexc_7PI, Colors[pulse_value]+'--', label=r'$\sigma_7=$'+str('%.4g'%SevenPhotonAbsorption_params[0][0])+tau_text) #NOTE adds value of sigma_3 here
        del Nexc_7PI, SevenPhotonAbsorption_params, Intensity_sorted_mask
        return plot_field
    #}}}
    
    #NOTE: the photon number does NOT change monotoneously with increasing intensity
    #TODO: one has to code a fitting method that sweeps intensity, but accepts ANY number of photons on each interval. 
    #NOTE: Maybe use the computed error to accept or reject each fitting? One can try all numbers of photons on each interval. 
    
    if(pulseWavelength==800e-9): 
        SevenPhotonIntensityRegime_min = 0e0
        SevenPhotonIntensityRegime_max  = 0E0
        
        SixPhotonIntensityRegime_min = 0e0
        SixPhotonIntensityRegime_max  = 0E0
        
        FivePhotonIntensityRegime_min = 0e0
        FivePhotonIntensityRegime_max  = 0E0
        
        FourPhotonIntensityRegime_min = 0e0
        FourPhotonIntensityRegime_max = 0e0 #W/m2
        
        ThreePhotonIntensityRegime_min = 0e0
        ThreePhotonIntensityRegime_max = 0e0 #W/m2
        
        TwoPhotonIntensityRegime_min = 0e0
        TwoPhotonIntensityRegime_max = 1.4e15 #W/m2
        
        OnePhotonIntensityRegime_min = 1.4e15
        OnePhotonIntensityRegime_max = 1.4e17 #W/m2
        
        ListOfInterestingGammas = [0.10, 1.0, 10.]
        
        
    elif(pulseWavelength==1600e-9): 
        
        SevenPhotonIntensityRegime_min  = 0e0
        SevenPhotonIntensityRegime_max  = 0E0
        
        SixPhotonIntensityRegime_min  = 0e0
        SixPhotonIntensityRegime_max  = 0E0
        
        FivePhotonIntensityRegime_min  = 0e0
        FivePhotonIntensityRegime_max  = 0E0
        
        FourPhotonIntensityRegime_min = 1.1E14
        FourPhotonIntensityRegime_max = 1.6E15 #W/m2
        
        ThreePhotonIntensityRegime_min  = 0e0
        ThreePhotonIntensityRegime_max  = 0e0
        
        TwoPhotonIntensityRegime_min = 1.4E15
        TwoPhotonIntensityRegime_max = 3.2E16
        
        OnePhotonIntensityRegime_min = 3.2E16
        OnePhotonIntensityRegime_max = 1.6E17
        
        #TwoPhotonIntensityRegime_min = 3.4E15
        #TwoPhotonIntensityRegime_max = 4.7E16 #W/m2
        
        #OnePhotonIntensityRegime_min = 4.7E16
        #OnePhotonIntensityRegime_max = 1.3E17 #W/m2
        
        ListOfInterestingGammas = [0.10, 1.0, 10.]
        
    elif(pulseWavelength==3200e-9): 
        
        SevenPhotonIntensityRegime_min = 0E0 #2.1E14
        SevenPhotonIntensityRegime_max = 0E0 #1.3E15
        
        SixPhotonIntensityRegime_min = 0E0
        SixPhotonIntensityRegime_max = 0E0
        
        FivePhotonIntensityRegime_min = 2.1E14
        FivePhotonIntensityRegime_max = 1.3E15
        
        FourPhotonIntensityRegime_min = 0e0
        FourPhotonIntensityRegime_max  = 0E0
        
        ThreePhotonIntensityRegime_min = 1.6E15
        ThreePhotonIntensityRegime_max = 3.6E15 #W/m2
        
        TwoPhotonIntensityRegime_min = 3.61E15
        TwoPhotonIntensityRegime_max = 3.4E16 #W/m2
        
        OnePhotonIntensityRegime_min = 3.41E16
        OnePhotonIntensityRegime_max = 7.3E16 #W/m2
        
        ListOfInterestingGammas = [0.10, 1.0, 10.]
    else: 
        SevenPhotonIntensityRegime_min = 0e0
        SevenPhotonIntensityRegime_max  = 0E0
        
        SixPhotonIntensityRegime_min = 0e0
        SixPhotonIntensityRegime_max  = 0E0
        
        FivePhotonIntensityRegime_min = 0e0
        FivePhotonIntensityRegime_max  = 0E0 
        
        FourPhotonIntensityRegime_min = 0e0
        FourPhotonIntensityRegime_max  = 1.1E14
        
        ThreePhotonIntensityRegime_min = 0e0
        ThreePhotonIntensityRegime_max = 0e0 #W/m2
        
        TwoPhotonIntensityRegime_min = 0e0
        TwoPhotonIntensityRegime_max = 1.3e15 #W/m2
        
        OnePhotonIntensityRegime_min = 1e12
        OnePhotonIntensityRegime_max = 1.3e15 #W/m2
        
        ListOfInterestingGammas = [0.10, 1.0, 10.]
        
    omega = 2.*np.pi*c/wavelengths1[0]
    NumberOfPhotons = np.divide(Egap_SI/hbar, omega)
    print NumberOfPhotons
    NumberOfPhotons_round = np.ceil(NumberOfPhotons)
    print Header+"** INFO: NUMBER OF PHOTONS = "+str(NumberOfPhotons)+", rounded to: "+str(NumberOfPhotons_round)
    
    # Here we attempt to fit the excitation coefficient in the multiphotonic regime. 
    # Hence, it depends on the number of photons required for overcoming the direct band gap energy. 
    # NOTE: we employ only ONE band gap energy here, despite there is a full band structure to consider. 
    # Note that tunneling effect only depends on the intensity. Therefore it could look like
    # 1-photon-absorption. 
    
    try:
        #NOTE: assumption here is the MAIN BAND GAP is employed for the transition. But it may NOT be the case. 
        #plotMe = (pulse_value == 20E-15)
        #if(NumberOfPhotons_round >= 4 and NumberOfPhotons_round < 5): 
        print Header+"** Info: fitting with 7-photon absortion..."
        plot_field = FitSevenPhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, 7, SevenPhotonIntensityRegime_min, SevenPhotonIntensityRegime_max)
        plot1 += [plot_field]
    except: 
        print Header+"** Fit using 7-photon absorption failed. "
    
    try:
        #NOTE: assumption here is the MAIN BAND GAP is employed for the transition. But it may NOT be the case. 
        #plotMe = (pulse_value == 20E-15)
        #if(NumberOfPhotons_round >= 4 and NumberOfPhotons_round < 5): 
        print Header+"** Info: fitting with 6-photon absortion..."
        plot_field = FitSixPhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, 6, SixPhotonIntensityRegime_min, SixPhotonIntensityRegime_max)
        plot1 += [plot_field]
    except: 
        print Header+"** Fit using 6-photon absorption failed. "
    
    try:
        #NOTE: assumption here is the MAIN BAND GAP is employed for the transition. But it may NOT be the case. 
        #plotMe = (pulse_value == 20E-15)
        #if(NumberOfPhotons_round >= 4 and NumberOfPhotons_round < 5): 
        print Header+"** Info: fitting with 5-photon absortion..."
        plot_field = FitFivePhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, 5, FivePhotonIntensityRegime_min, FivePhotonIntensityRegime_max)
        plot1 += [plot_field]
    except: 
        print Header+"** Fit using 5-photon absorption failed. "
        
    try:
        #NOTE: assumption here is the MAIN BAND GAP is employed for the transition. But it may NOT be the case. 
        #plotMe = (pulse_value == 20E-15)
        #if(NumberOfPhotons_round >= 4 and NumberOfPhotons_round < 5): 
        print Header+"** Info: fitting with 4-photon absortion..."
        plot_field = FitFourPhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, 4, FourPhotonIntensityRegime_min, FourPhotonIntensityRegime_max)
        plot1 += [plot_field]
    except: 
        print Header+"** Fit using 4-photon absorption failed. "
    try:
    #elif(NumberOfPhotons_round >= 3 and NumberOfPhotons_round < 4):
        print Header+"** Info: fitting with 3-photon absortion..."    
        plot_field = FitThreePhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, 3, ThreePhotonIntensityRegime_min, ThreePhotonIntensityRegime_max)
        plot1 += [plot_field]
    except: 
        print Header+"** Fit using 3-photon absorption failed. "
    try:
    #elif(NumberOfPhotons_round >= 2 and NumberOfPhotons_round < 3): 
        print Header+"** Info: fitting with 2-photon absortion..."
        plot_field = FitTwoPhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, 2, TwoPhotonIntensityRegime_min, TwoPhotonIntensityRegime_max)
        plot1 += [plot_field]
    except: 
        print Header+"** Fit using 2-photon absorption failed. "
    try:
    #elif (NumberOfPhotons_round >= 1 and NumberOfPhotons_round < 2): 
        print Header+"** Info: fitting with 1-photon absortion..."
        plot_field = FitOnePhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, 1, OnePhotonIntensityRegime_min, OnePhotonIntensityRegime_max)
        plot1 += [plot_field]
    except: 
        print Header+"** Fit using 1-photon absorption failed. "
        #else: 
            #print Header+"** Warning: could not fit anything, surely because you are using the 2-band vision! :) "
            
    #except: 
        #print Header+"** Fit using N-photon absorption failed. "
    # NOTE: this is NOT 3-photon absorption regime, this is NON-PERTURBATIVE regime. Transition to tunneling. 
    
    print Header+"** Info: fitting the non-perturbative regime..."
    TunnelingIntensityRegime_min = 2E15 #TwoPhotonIntensityRegime_max #1e12
    TunnelingIntensityRegime_max = 2E17  #1.3e15 #W/m2
    #try:
        #print Header+"** Info: fitting the tunneling rate..."
        #plot_field = FitOnePhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, 1, TunnelingIntensityRegime_min, TunnelingIntensityRegime_max)
        #plot1 += [plot_field]
    #except: 
        #print Header+"** Fit using tunneling rate failed. "
    #try: 
        #plot_field = FitTwoPhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, 2, ThreePhotonIntensityRegime_min, ThreePhotonIntensityRegime_max)
        #plot1 += [plot_field]
    #except: 
        #print Header+"** Warning: 3-ph fitting of the data failed for tau="+str(pulse_value*1E15)+" fs."
    for gammaValues in ListOfInterestingGammas:
      ax1.axvline(x=np.multiply(xscale, IntensityAtGamma(gammaValues, Egap_SI, meff_Waroquiers_LDA, pulseWavelength)), color='k', linestyle='--', linewidth=0.5)
    ax1.axhline(y=  Density, color='k', linestyle='--', linewidth=0.25)
    #ax1.axhline(y=2*Density, color='k', linestyle='--', linewidth=0.25)
    #ax1.axhline(y=3*Density, color='k', linestyle='--', linewidth=0.25)
    ax1.axhline(y=4*Density, color='k', linestyle='--', linewidth=0.25)
  #}}}
  ax1.set_xlabel(r"Intensity (W/cm$^{2}$)")
  ax1.set_ylabel(r"Final exc. density (m$^{-3}$)")
  #pplt.xlabel(r'Peak intensity (W/m$^2$)')
  #pplt.ylabel(r'Final exc. density (m$^{-3}$)')
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  
  labels1 = [l.get_label() for l in plot1]
  ax1.legend(plot1, labels1, loc='best', fancybox=True, framealpha=1, fontsize=OneSizeFitsAll)
  
  pplt.tight_layout()
  #print Header+"Adding the adjusted coefficients to check our fitting law."
  
  ax2 = pplt.subplot(plotB)
  if(MoveIon): 
      MoveIonTitle=", w/ phonons."
      MoveIonFile="MoveIon"
  else: 
      MoveIonTitle=""
      MoveIonFile=""
  pplt.title(r'$\lambda = $'+str(int(pulseWavelength*1E9))+r' nm, $n=\left\lceil\frac{E_g^\Gamma}{\hbar \omega}\right\rceil=$'+str(int(NumberOfPhotons_round))+MoveIonTitle)
  print Header+"Generating comparison with Keldysh-Gruzdev model..." ## Looking forward to plot w_PI[intensity] from Keldysh_Gruzdev theory. 
  order = 50
  Efieldmin = np.min(Efield1_uniq)
  Efieldmax = 5E10 # np.max(Efield1_uniq)
  Efield1_Keldysh_log = np.linspace(np.log10(Efieldmin), np.log10(Efieldmax), 2000)
  Efield1_Keldysh     = np.power(10., Efield1_Keldysh_log)
  print Efield1_Keldysh
  
  for gammaValues in ListOfInterestingGammas:
      ax2.axvline(x=xscale*IntensityAtGamma(gammaValues, Egap_SI, meff_Waroquiers_LDA, pulseWavelength), color='k', linestyle='--', linewidth=0.5)
#  pplt.ylim((1E16,4.*5e28))
  #ax1.grid()
  #if(PlotWpiAxis==True): #{{{
  plot2=[]
  for pulse_value in pulses_uniq: #{{{
      # select the corresponding wavelength in results
      PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses, pulse_value, IndexTable['pulses1'])
      PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses_pulses, pulse_value, IndexTable['pulses2']) #keeps only pulses of the same duration.

      runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_pulses)
      
      lists = sorted(itertools.izip(*[FieldToIntensity(Efield1), Ne_final, pulses1]))
      Intensity_sorted, Ne_final_sorted, pulses1_sorted = list(itertools.izip(*lists))
      
      plot21, = ax2.plot(np.multiply(xscale, Intensity_sorted), np.divide(Ne_final_sorted, np.multiply(1.,pulses1_sorted)), Colors[pulse_value]+Symbols[pulse_value], label=r'$\tilde{w}_{\mathrm{PI}}^{\mathrm{TDDFT}}$, $\tau=$'+str(int(pulse_value*1E15))+' fs')
      plot2 += [plot21]
  #}}}
  
  if(pulseWavelength==800e-9): 
    EffectiveMassList=[meff_Waroquiers_LDA] #[0.05]
    #EffectiveMassList=[meff_Waroquiers_LDA, 1.] #, meff_Waroquiers_LDA]
  elif(pulseWavelength==1600e-9):
    EffectiveMassList=[meff_Waroquiers_LDA]
  elif(pulseWavelength==3200e-9):
    EffectiveMassList=[meff_Waroquiers_LDA] #[0.025] #, meff_Waroquiers_LDA]
  else: 
    EffectiveMassList=[0.01, 0.025, 0.05, 0.1, meff_Waroquiers_LDA]
  for meff_Keldysh in EffectiveMassList: #{{{
    PeakWpi_Gruzdev = GenerateKeldyshDatabase(Egap_SI, meff_Keldysh, pulseWavelength, Efield1_Keldysh*np.sqrt(OpticalIndex), order, 1.) #NOTE: optical index should be treated as such, not through adiabadicity coefficient, or? 
    
    Acoeff=1
    PeakWpi_Atomic_Keldysh = IonizationRateAtoms(pulseWavelength, meff_Keldysh, Egap_SI, Efield1_Keldysh, order, 1., Ntotal) #Optical index in gas in 1. 
    print Header+"Info: Atomic Keldysh: "+str(np.min(PeakWpi_Atomic_Keldysh))+", "+str(np.max(PeakWpi_Atomic_Keldysh))
    #PeakWpi_Gulley = GenerateKeldyshGulleyDatabase(Egap_SI, meff_Keldysh, pulseWavelength, Efield1_Keldysh, order, 1.)
    gamma = gammaKeldysh(Egap_SI, meff_Keldysh, Efield1, wavelengths1, 1.) #for gamma, ONE MUST NOT USE OPTICAL INDEX (because gamma is oscillation time divided by tunneling time, hence optical index must vanish)
    #ax2x2 = ax2.twiny()
    # We need to plot values of intensities for which gamma=0.1, 1, 3, 5, 10.
    # d e f IntensityAtGamma(gamma, Egap_SI, meff, wavelength):
        #omega = 2.*np.pi * c / wavelength
        #Ipeak = meff*m_e*Egap_SI*omega**2*c*epsilon_0 / (2. * e**2 * gamma ** 2)
        #return Ipeak
        
    ## GRUZDEV    
    ## NOTE: Gruzdev treak consists in plotting (external_field, w_pi[internal_field]).      
    plot22, = ax2.plot(np.multiply(xscale, FieldToIntensity(Efield1_Keldysh)), PeakWpi_Gruzdev, '-k', label=r'$w_{\mathrm{PI}}^{\mathrm{KG}}$') #, $m_{\mathrm{eff}}=$'+str(np.round(meff_Keldysh,3))+'.')
    plot2    += [plot22]
    
    # ATOMIC KELDYSH
    print Header+"Preparing atomic Keldysh..."
    plot23, = ax2.plot(np.multiply(xscale, FieldToIntensity(Efield1_Keldysh)), PeakWpi_Atomic_Keldysh, '-.k', label=r'$w_{\mathrm{PI}}^{\mathrm{at-anal}}$') #, $m_{\mathrm{eff}}=$'+str(np.round(meff_Keldysh,3))+'.')
    plot2    += [plot23]
  #}}}  
  ## ZHUKOV
  PathPrefix="/home/hilase/Documents/spp-extended-theory/Keldysh/"
  PeakWpi_Zhukov = VZ_generateWpiTables(Efield1_Keldysh, 0., pulseWavelength, 800e-9, 0., 0., Egap_SI, meff_Keldysh, Ntotal, PathPrefix) #we neglect optical index in atomic theory
  
  try: #we plot as function of external field
      plot23, = ax2.plot(np.multiply(xscale, FieldToIntensity(Efield1_Keldysh[np.nonzero(PeakWpi_Zhukov)])), PeakWpi_Zhukov[np.nonzero(PeakWpi_Zhukov)], ':k', label=r'$w_{\mathrm{PI}}^{\mathrm{at-num}}$') #, $m_{\mathrm{eff}}=$'+str(np.round(meff_Waroquiers_LDA,3))+'.') #X axis: external intensity. #Y axis: internal (compressed) intensity. 
      plot2    += [plot23]
  except: 
      print Header+"** Warning: no data available from VP Zhukov for this wavelength."
      
    ## GULLEY 
    #plot23, = ax2.loglog(FieldToIntensity(Efield1_Keldysh), PeakWpi_Gulley, '-.', label=r'$w_{\mathrm{PI}}^{\mathrm{Gu.}}$, $m_{\mathrm{eff}}=$'+str(np.round(meff_Keldysh,3))+'.')
    #plot2    += [plot23]
    
  ## TUNNELING SECTION
  
  Egap_SI_Tunnel = Egap_SI #direct band gap for LDA. Indirect requires phonons (or local field effects)
  
  EnableStarkEffect=False
  if(pulseWavelength  == 800e-9): 
    EffectiveMassList=[meff_Waroquiers_LDA]
  elif(pulseWavelength==1600e-9):
    EffectiveMassList=[meff_Waroquiers_LDA]
  elif(pulseWavelength==3200e-9):
    EffectiveMassList=[meff_Waroquiers_LDA]
  else: 
    EffectiveMassList=[meff_Waroquiers_LDA]
  wPI_min={800e-9: 1E36, 1600e-9: 1E35, 2200e-9: 1E30, 2600e-9: 1E30, 3200e-9: 1E34}
  
  for meff_Tunnel in EffectiveMassList: #meff_Waroquiers_LDA, meff_VeryHighField]: 
  #{{{ #NOTE Stark effect is enabled here, not for the wavelength approach. 
    PeakWpi_Tunneling   = np.clip(KeldyshTunnelingLimit(Egap_SI_Tunnel, meff_Tunnel, pulseWavelength, Efield1_Keldysh*np.sqrt(OpticalIndex), EnableStarkEffect, 1.), 1E30,2E46)
    #NOTE: tunneling barrier does not change as function of refractive index. Hence, we remove the optical index in case of tunneling limit. 
    plot23, = ax2.plot(np.multiply(xscale, FieldToIntensity(Efield1_Keldysh)), PeakWpi_Tunneling[np.nonzero(PeakWpi_Tunneling)], '--k', label=r'$w_{\mathrm{PI}}^{\mathrm{Tunnel}}$') #, $m_{\mathrm{eff}}=$'+str(np.round(meff_Tunnel,3))+'.')
    plot2    += [plot23]
    
    ## TUNNELING WITH INDIRECT BAND GAP.
    #PeakWpi_Tunneling_i   = np.clip(KeldyshTunnelingLimit(Egap_SI_indirect, meff_Tunnel, pulseWavelength, Efield1_Keldysh, EnableStarkEffect, 1.), 1E30,1E46)
    #plot24, = ax2.loglog(FieldToIntensity(Efield1_Keldysh), PeakWpi_Tunneling_i, label=r'$w_{\mathrm{TI}}^{i}$, $m_{\mathrm{eff}}=$'+str(np.round(meff_Tunnel,3))+'.')
    #plot2    += [plot24]
  #}}}
  pplt.ylim((wPI_min[pulseWavelength],1E45)) #we adjust min of the wPI as function of wavelength
  # Interpolation of a joint law of photoionization
  #PeakWpi_Effective=np.power(np.power(PeakWpi_Tunneling,-1)+np.power(PeakWpi_Gruzdev,-1),-1)
  
  #plot23, = ax2.loglog(FieldToIntensity(Efield1_Keldysh), PeakWpi_Effective, 'k-', linewidth=2, label=r'$w_{\mathrm{eff}}$')
  #plot2    += [plot23]
  
  labels2 = [l.get_label() for l in plot2]
  ax2.legend(plot2, labels2, loc='best', fancybox=True, framealpha=1, fontsize=OneSizeFitsAll)
  #ax12.set_yscale('log')
  #ax12.legend(loc='best', fancybox=True, framealpha=0)
  ax2.set_xlabel(r'Intensity (W/cm$^2$)')
  ax2.set_ylabel('$w_{\mathrm{PI}}$ (m$^{-3}$.s$^{-1}$)')
  ax2.set_xscale('log') #compared to using loglog(), the xtics and ytics are better for publication. 
  ax2.set_yscale('log')
  pplt.tight_layout()
  #ax2.legend(plot1, labels1, loc='best', fancybox=True, framealpha=0)
  #ax2.axes.get_xaxis().set_visible(False)
  #ax2.axes.get_yaxis().set_visible(False)
  #else: 
  
  ## This should be the inset that we would put in Fig. 3a. 
  # Previously, it was prepared using GNUplot file plotFieldWpiPeaks.gnu.
  
  #if(PlotVirtualExcitation): #{{{ #NOTE: no longer used, since it has no meaning. Gauge dependent quantity. 
    
    #ax3 = pplt.subplot(plotC)
    #plot1 = [] #np.empty((0,0))
    #for pulse_value in pulses_uniq: #{{{  
        ##pplt.title('Virtual exc. factor')
        ## select the corresponding wavelength in results
        #PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses, pulse_value, IndexTable['pulses1'])
        #PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses_pulses, pulse_value, IndexTable['pulses2']) #keeps only pulses of the same duration.

        #runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, ElectronEnergy_max, Energy_max, W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_pulses)
        
        #print Header+"** Info: saving the data in CSV file for fitting. "
        #filename = "EffectOfPulseIntensity-Duration-"+str(pulseWavelength*1E9)+"nm-tau"+str(int(pulse_value*1E15))+"fs"
        #print Header+"** Info: saving in "+filename+".csv"
        ##ExportToTxt(PlottingSinglePulses_pulses, filename+".csv")
        #lists = sorted(itertools.izip(*[FieldToIntensity(Efield1), irreversibility, pulses1]))
        #Intensity_sorted, irreversibility_sorted, pulses1_sorted = list(itertools.izip(*lists))
        
        #plot_field, = ax2.plot(1E-4*Intensity_sorted, irreversibility_sorted, Colors[pulse_value]+'-o', label=r'$\zeta$, $\tau=$'+str(round(pulse_value*1E15))+' fs')
        #plot1 += [plot_field]
    ##}}}
    #pplt.xlabel(r'Peak intensity (W/cm$^2$)')
    #pplt.ylabel(r'Virtual exc. $\zeta$')
    #pplt.ylim((0,1))
    #ax3.set_xscale('log')
    #ax3.grid()
    
    #ax32 = ax3.twinx()
    #for pulse_value in pulses_uniq:
      ## select the corresponding wavelength in results
      #PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses, pulse_value, IndexTable['pulses1']) #keeps only 800 nm wavelenghts.

      #runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, ElectronEnergy_max, Energy_max, W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_pulses)
      
      #plot_field, = ax32.plot(FieldToIntensity(Efield1), gammaKeldysh(Egap_SI, meff_Waroquiers_LDA, Efield1, wavelengths1), 'k-', label=r'$\gamma$')
      ##plot1 += [plot_field]
    #pplt.ylabel(r'Adiabadicity $\gamma$')
    ##pplt.ylim((0,1))
    #ax32.set_yscale('log')
      
    ## ax2.set_yscale("log")
    
    #labels3 = [l.get_label() for l in plot3]
    ##ax22.legend(plot1, labels1, loc='center right', fancybox=True, framealpha=0)
    ##ax1.legend(loc=1, borderaxespad=0.)
    ## ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ##ax22.legend(loc=2, borderaxespad=0.)
    #pplt.tight_layout()
  ##}}}
  
  
  def FitDrudeModel(n_exc, nu):
    #omega       = 2.*np.pi*c/wavelength
    omega_p2    = np.multiply(n_exc, e**2 / m_e / meff / epsilon_0)
    Epsilon = 1.  - np.divide(omega_p2, np.power(omega, 2)) * np.divide(1.,1. + np.divide(1.j * nu, omega)) #conduction band absorption, no interband contribution
    return Epsilon

  # d e f FitDrudeModelReal(n_exc, nu):
      #Epsilon = FitDrudeModel(n_exc, nu)
          ##1.  - np.divide(omega_p2, np.power(omega, 2)) * np.divide(1.,1. + np.divide(1.j * nu, omega)) #conduction band absorption, no interband contribution
      ## the coupling/interference between interband and intraband motions is accounted for in the complete dielectric function!
      #DrudeEnergy          = SuperCellVolume_SI / e * epsilon_0 * omega * Epsilon.imag * np.multiply(Efield_X_SI_interp_Ne, Efield_X_SI_interp_Ne)
      #TotalDrudeEnergy     = np.multiply(np.cumsum(DrudeEnergy), dt_deriv[0])+ExcessEnergy
      
      #return TotalDrudeEnergy

  
  if(PlotEnergy):
    ax3 = pplt.subplot(plotC)
    #pplt.title('Virtual exc. factor')
    
    plot3 = []
    for pulse_value in pulses_uniq: #{{{
        # select the corresponding wavelength in results
        PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses, pulse_value, IndexTable['pulses1'])
        PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses_pulses, pulse_value, IndexTable['pulses2']) #keeps only pulses of the same duration.

        runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_pulses)
        
        lists = sorted(itertools.izip(*[FieldToIntensity(Efield1), ElectronEnergy_max, pulses1, Ne_final]))
        Intensity_sorted, Energy_final_sorted, pulses1_sorted, Ne_final_sorted = list(itertools.izip(*lists))
        #print Header+"Ready to plot: "
        
        plot_field, = ax3.plot(np.multiply(xscale, Intensity_sorted), Energy_final_sorted, Colors[pulse_value]+'-'+Symbols[pulse_value], label=r'$\tau=$'+str(int(pulse_value*1E15))+" fs")
        plot3 += [plot_field]
        ## TODO: can we adjust collision frequency on Drude here? 
        
        #try: 
        #FitDrudeModel_parameters = curve_fit(FitDrudeModelReal, Ne_final_sorted, Energy_final_sorted, bounds=(1E13,1E16))
        # We have to minimize ElectronEnergy - 
        #print Header+"FITTING DRUDE MODEL."
        #print np.divide(1.,FitDrudeModel_parameters)
        #nu_t = FitDrudeModel_parameters[0]
        
        #TODO: still have to use epsilon_Drude to compute Energy. 
        #plot_field, = ax3.plot(Intensity_sorted, FitDrudeModel(Ne_final_sorted,nu_t), Colors[pulse_value]+'-.', label=r'Drude, $\tau=$'+str(int(pulse_value*1E15))+r" fs, $\nu^{-1}=$"+str(round(nu_t,2)))
        #plot3 += [plot_field]
        #except: 
        #print Header+"** WARNING: auto-fit of collision frequency failed."
        #nu_t = []
        
        
    #}}}
    labels3 = [l.get_label() for l in plot3]
    ax3.legend(plot3, labels3, loc='best', fancybox=True, framealpha=1, fontsize=OneSizeFitsAll)
    pplt.xlabel(r'Peak intensity (W/cm$^2$)')
    pplt.ylabel(r'Energy $\xi_{\mathrm{el}}^{\mathrm{TDDFT}}$ (eV)')
    #pplt.ylim((0,1))
    #locmin = mp.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1, 2, 3,4,5,6,7,8,9,10)) 
    #locmin = mp.ticker.LogLocator(base=10.0,numticks=4)
    #ax1.get_xaxis().set_major_locator(mp.ticker.LogLocator(base=10.0,numticks=4))
    #ax1.get_xaxis().set_minor_locator(mp.ticker.LogLocator(base=10.0,numticks=50))
    #ax1.xaxis.set_minor_locator(mp.ticker.AutoMinorLocator())
    #ax1.get_xaxis().set_minor_formatter(mp.ticker.NullFormatter())
    #ax1.get_xaxis().set_minor_formatter(mp.ticker.ScalarFormatter())
    
    
    #ax2.xaxis.set_minor_locator(locmin)
    #ax2.xaxis.set_minor_formatter(mp.ticker.NullFormatter())
    
    #ax3.xaxis.set_minor_locator(locmin)
    #ax3.xaxis.set_minor_formatter(mp.ticker.NullFormatter())
    #ax3.yaxis.set_major_locator(locmin)
    #ax3.yaxis.set_minor_formatter(mp.ticker.NullFormatter())
    
    ax3.set_xscale('log')
    #ax3.set_yscale('symlog', linthreshy=1e-5)
    ax3.set_yscale('log')
    
    pplt.tight_layout()
    #ax3.grid()
    #ax32 = ax3.twinx()
    #for pulse_value in pulses_uniq:
      # select the corresponding wavelength in results
      #PlottingSinglePulses_pulses = FilterDatabase(PlottingSinglePulses, pulse_value, IndexTable['pulses1']) #keeps only 800 nm wavelenghts.

      #runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, ElectronEnergy_max, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses_pulses)
      
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
  ExportName = "EffectOfPulseIntensity-Duration-"+str(pulseWavelength*1E9)+"nm"+MoveIonFile
  pplt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.eps')
  pplt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.png')
  print Header+"** Saved: 'Metaplot-SinglePulse-SingleColor-'"+ExportName+".png."
  if(ShowFinalPlot==1):
    pplt.show()

  del PlottingSinglePulses, wavelength_uniq, CEP1_uniq
  del ax2
  if(PlotEnergy): 
      del ax3
#}}}

## Single-pulse: Final electron density as function of intensity. Each curve is for a different CEP. 
# TODO: current and energy, and harmonic spectrum as function of CEP.
def plotEffectOfIntensity_CEP(PlottingTable, pulseWavelength=800e-9, pulseDuration=20e-15, FilterCEP=False, CEP0=0., PulseShape=3, pulseRiseTime=0, ShowFinalPlot=0, Egap_SI=2.56*e): #{{{
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
  pplt.title(r'$n_{\mathrm{exc}}(I)$, $\tau=$'+str(int(pulseDuration*1E15))+' fs.')
  # Number of sweeped parameters: 3 (Efield, Delay, CrossPolarization)
  #pplt.scatter(FieldToIntensity(Efield1), Ne_final, c=Delay, s=(30.+CrossPolarization*50.), label=r'$N_{exc}$[Field]')
  plot1 = [] #np.empty((0,0))
  for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
    PlottingSinglePulses_CEP = FilterDatabase(PlottingSinglePulses, CEP1_value, IndexTable['CEP1']) #select the local CEP
    
    ## prints the current state of selection and stores the values of sweeped parameters
    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_CEP)
    
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
    print Header+"{Fore.RED}** Debug info: "+str(np.shape(plot_field))
    plot1 += [plot_field] # does cumulate each CEP now. Type: Line2D, not np.array. 
    
    #ax1.legend(loc=4, borderaxespad=0.) #bbox_to_anchor=(1.05, 1), 
  labels1 = [l.get_label() for l in plot1]
  ax1.legend(plot1, labels1, loc=2)
  plt.xlabel(r'Peak intensity (W/cm$^2$)')
  plt.ylabel(r'Final density (m$^{-3}$)')
  #  pplt.ylim((1E16,4.*5e28))
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  #plot1_list = plot1.tolist()
    #print Efield1, CEP1_value
    
  #Adding axis of gammaKeldysh, but for all CEPs !!
  #runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, ElectronEnergy_max, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses)
  
  ax12 = ax1.twinx()
  for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
    PlottingSinglePulses_CEP = FilterDatabase(PlottingSinglePulses, CEP1_value, IndexTable['CEP1']) #select the local CEP
    
    ## prints the current state of selection and stores the values of sweeped parameters
    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_CEP)
    
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

## Nexc ( CEP )  
def plotEffectOfCEP(PlottingTable, pulseWavelength=800e-9, pulseDuration=20e-15, FieldAmplitude=1E9, PulseShape=3, pulseRiseTime=0, ShowFinalPlot=0, Egap_SI=2.56*e): #{{{
  print "#====== Plot: SINGLE PULSE, EFFECT OF THE CEP, SAME WAVELENGTHS, SAME POLARIZATIONS: final density as function of field ======"
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
  pplt.title(r'$n_{\mathrm{exc}}(I)$, $\tau=$'+str(int(pulseDuration*1E15))+r' fs, $E_{peak}=$'+str(round(FieldAmplitude*1E-9))+' V/nm.')
  # Number of sweeped parameters: 3 (Efield, Delay, CrossPolarization)
  #pplt.scatter(FieldToIntensity(Efield1), Ne_final, c=Delay, s=(30.+CrossPolarization*50.), label=r'$N_{exc}$[Field]')
  plot1 = [] #np.empty((0,0))
  #for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
  #PlottingSinglePulses_CEP = FilterDatabase(PlottingSinglePulses, CEP1_value, IndexTable['CEP1']) #select the local CEP
  
  ## prints the current state of selection and stores the values of sweeped parameters
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses)
  
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
  print Header+"{Fore.RED}** Debug info: "+str(np.shape(plot_field))
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
  #runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, ElectronEnergy_max, Energy_max, W_pi, PulseType, RiseTime = ExtractData(PlottingSinglePulses)
  
  ax12 = ax1.twinx()
  #for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
  #PlottingSinglePulses_CEP = FilterDatabase(PlottingSinglePulses, CEP1_value, IndexTable['CEP1']) #select the local CEP
  
  ## prints the current state of selection and stores the values of sweeped parameters
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses)
  
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

## (Final density | Virtual excitation rate) (Intensity) with comparing FreezeHXC if available.
def plotEffectOfFrozenXC(PlottingTable): #{{{
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
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_CEP_NoFreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), Ne_final, 'ow', label=r'$N_{exc}^{final}$') #, CEP='+str(round(CEP1_value*180./pi)))
  
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_CEP_FreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), Ne_final, 'ok', label=r'$N_{exc}^{final}$, FrozenHXC')#, CEP='+str(round(CEP1_value*180./pi)))
  
  ax1.legend(loc='best', fancybox=True, framealpha=0) #bbox_to_anchor=(1.05, 1), 
  #pplt.xlabel(r'Peak intensity (W/m$^2$)')
  pplt.ylabel(r'Final density (m$^{-3}$)')
  #  pplt.ylim((1E16,4.*5e28))
  ax1.set_xscale('log')
  ax1.set_yscale('log')

  #Adding axis of gammaKeldysh
  ax12 = ax1.twinx()
  
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_CEP_NoFreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), Ne_final/pulses1, 'ow') #, label=r'$N_{exc}/\tau$') #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
  
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_CEP_FreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), Ne_final/pulses1, 'ok') #, label=r'$N_{exc}/\tau$, FrozenHXC') #, I='+str(FieldToIntensity(Efield1_value))+r'
  
  ax12.set_yscale('log')
  ax12.set_ylabel(r'$w_{PI}$ (m$^{-3}$.s$^{-1}$)')

  ax2 = pplt.subplot(plotB)
  pplt.title(r'Virtual excitation $\zeta$')
  
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_CEP_NoFreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), irreversibility, 'ow', label='Virtual exc.')
  
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_CEP_FreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), irreversibility, 'ok', label='Virtual exc., FreezeHXC')
  
  pplt.xlabel(r'Peak intensity (W/cm$^2$)')
  pplt.ylabel(r'Virtual exc. $\zeta$')
  pplt.ylim((0,1))
  ax2.set_xscale('log')
  #adding the Keldysh adiabatic factor
  ax22 = ax2.twinx()
  
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_CEP_NoFreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), gammaKeldysh(Egap_SI, meff, Efield1, wavelengths1), 'or', label=r'Adiabaticity $\gamma$') #BUG: 2.56 must not be specified by hand here! 
  
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_CEP_FreezeHXC)
  
  pplt.plot(FieldToIntensity(Efield1), gammaKeldysh(Egap_SI, meff, Efield1, wavelengths1), 'or') #, label='Adiabaticity') #BUG: 2.56 must not be specified by hand here! 
  
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

## Plots (Nexc | Wpi | Energy) (wavelength) (1D plot) by introducing <FieldToTolerance> percents tolerance on the electric field. 
# This enables to analyze the role of laser wavelength on the final excited electron density, excitation rate and laser-induced energy increase. Provides predictions on the Stark effect. 
# Wpi is compared to analytical Keldysh theory. 
# Energy is compared with analytical Drude model. 
def plotEffectOfWavelength_Tolerance(PlottingTable, FieldTolerance=0.05, pulseDuration=20e-15, CEP=0., PulseType=3, RiseTime=2, ShowFinalPlot=0, Egap_SI=2.56*e, Egap_SI_indirect=0.51*e, wavelengths1_min=200e-9, wavelengths1_max=3300e-9, ShowDensities=False, PlotIrreversibility=False, PlotFreeCarrierAbsorption=True, EnergyAxis=True): #{{{
  PlotStarkShift=False
  PlotValleysInset=False #True: plot valleys, false: plot hills (max of Wpi)
  AccountForOpticalIndex=False #trying without optical index
  global n0
  
  meff_Waroquiers_LDA = 0.2226 #effective mass of electron
  #[Laflamme Janssen, J.; Gillet, Y.; Poncé, S.; Martin, A.; Torrent, M. & Gonze, X. Precise effective masses from density functional perturbation theory Physical Review B, 2016, 93, 205147]
  #meff = 0.14 #0.14: fits for Si in MPI regime at 800 nm.
  #Field dependent fitting of effective mass
  meff_LowField = meff_Waroquiers_LDA #0.14
  meff_HighField = meff_Waroquiers_LDA #0.14: fits for Si in MPI regime at 800 nm.
  meff_VeryHighField = meff_Waroquiers_LDA #0.5
  
  meff_Tunnel = meff_Waroquiers_LDA #fits well if using INDIRECT band gap for 800 nm
  print "#====== Plot: SINGLE PULSE, EFFECT OF THE WAVELENGTH, SAME CEP, SAME POLARIZATIONS: final density as function of field ======"
  #Filter: double pulses, same wavelengths, same polarization
  
  LineWidthCoeff_abs = 0.5 #controls absolute thickness of all curves
  LineWidthCoeff_rel = 1.0 #controls thickness of Keldysh-Gruzdev curves relatively to TDDFT
  antialias=True
  
  #FieldOfInterest = [5e8, 1e9, 2e9, 3E9, 4E9, 5e9, 1e10]
  FieldOfInterest = [5e8, 7.5E8, 1e9, 1.5e9, 2e9, 3e9, 4E9, 5e9] #, 1e10]
  NumberOfColors=len(FieldOfInterest)
  ColorList=[plt.cm.seismic(i) for i in np.linspace(0, 1, NumberOfColors)] 
  # Not bad! RdYlBu, coolwarm, seismic
  # Too soft: Copper, Spectral
  # Not perfect. RdBu. RdYlGn.
  # Too dynamic: bwr
  # To test: , brg, rainbow
  #Colors={5e8: 'c', 7.5e8: 'r', 1e9: 'b', 1.5e9: 'g', 2e9: 'y', 3e9: 'm'} #, 1e10: 'k'}
  Colors={FieldOfInterest[i]:ColorList[i] for i in range(0,NumberOfColors)}
  LineWidths={FieldOfInterest[i]:np.arange(1.5,4,NumberOfColors) for i in range(0,NumberOfColors)}
  #EffectiveMass={1e9: meff_LowField, 2e9: meff_LowField, 5e9: meff_HighField, 1e10: meff_VeryHighField}
  EffectiveMass={FieldOfInterest[i]:meff_HighField for i in range(0,NumberOfColors)}
  # Simple pulse, single color
  
  def ImportPalikDatabase_epsilon_fromNK(wavelengths, DataFile="Si-Palik"): #{{{
    DataArray           = np.loadtxt(DataFile, delimiter="\t", skiprows=4)
    wavelengths_source  = DataArray[:,0]
    if(np.min(np.diff(wavelengths_source))==0e0): 
        print Header+"Error in wavelengths_source. Careful! Palik book combines different optical data, hence there can be overlap, sometimes! It must be removed to be imported properly here. " 
        #print np.diff(wavelengths_source)
        exit()
    OpticalN            = DataArray[:,1]
    OpticalK            = DataArray[:,2]
    order = 1
    fni = InterpolatedUnivariateSpline(wavelengths_source*1E-10, OpticalN, k=order)
    fki = InterpolatedUnivariateSpline(wavelengths_source*1E-10, OpticalK, k=order)
    OpticalN_new = fni(wavelengths)
    OpticalK_new = fki(wavelengths)
    epsilon = np.power(OpticalN_new+np.multiply(1j, OpticalK_new), 2)
    #epsR = np.real(epsilon); epsC = np.imag(epsilon)
    return epsilon
  #}}}

  ### TODO: MOVE TO libPlotting.py to access it for plotTimeKeldysh and plotTimeHHG.py. 
  def ImportPalikDatabase_epsilon_twofiles_eV(wavelengths, DataFile='Database/Si-TDLDA'): #{{{
    Header="** ImportPalikDatabase_epsilon_twofiles: "
    nfile = DataFile+"-epsR.csv" #TODO: rename nfile to ReEpsFile
    kfile = DataFile+"-epsC.csv" #TODO: rename kfile to ImEpsFile
    
    print nfile, kfile
    
    print "** Info: opening "+nfile+" and "+kfile+"."
    narray = loadtxt(nfile, delimiter="\t", skiprows=2)
    karray = loadtxt(kfile, delimiter="\t", skiprows=2)
    
    unit = 1E0 #Unit of the wavelength found in databases <nfile> and <kfile>. 
    
    # import wavelength, epsilonR, epsilonC
    energies1 = narray[:,0]*unit #/!\ its an energy!
    energies2 = karray[:,0]*unit #/!\ its an energy!
    n = narray[:,1]
    kk = karray[:,1]
    #numrows = 10000 #Number of rows to interpolate the data on. 
    base = 10
    # interpolate on new wavelength mesh using 1st order
    order=1 #linear interpolation
    #print Header+"** Info: Generating the new wavelength mesh: ("+str(np.amin(energies2))+", "+str(np.amax(energies2))+")"
    #wavelengths = np.logspace(np.amin(np.log10(energies2)), np.amax(np.log10(energies2)), num=numrows, base=base, endpoint = True)
    
    #print Header+"New wavelength mesh has "+str(numrows)+" rows."
    #print wavelengths
    
    print Header+"** Info: definition of interpolation functions..."
    energy1=np.sort(energies1)
    energy2=np.sort(energies2)
    del energies1, energies2 
    fni = InterpolatedUnivariateSpline(energy1, n, k=order)
    fki = InterpolatedUnivariateSpline(energy2, kk, k=order)
    
    ##Interpolation for one optical constant
    #try:
      #ni = fni(wavelength); ki = fki(wavelength)
      ##print wavelength, ni, ki
      #epsilon = (ni+1j*ki) #NOTE: we are picking up the epsRe, and epsIm directly here
      #print "" 
      #print "Interpolated permittivity at "+str(wavelength*1E9)+" nm = "+str(epsilon)
    #except: 
      #print "Interpolation for "+str(wavelength*1E9)+" nm failed."
	    
	    
    # We need to express the required energy mesh into a wavelength
    energies = hbar*np.divide(2.*np.pi * c, wavelengths) / e #eV
    
    ## WARNING: Names are misleading here: we actually work with dielectric permittivities! 
    # defining the new epsR and epsC on a common mesh
    ni = fni(energies) #eV
    ki = fki(energies)
    epsilons = np.add(ni, np.multiply(1.j, ki)) 
    return epsilons
  #}}}
  
  
  
  # TODO: Filtering function
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
  PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, wavelengths1_max, IndexTable['wavelengths1']) #we will try to add data from 3200 nm.

  ## TODO: check if CEP could count at smaller wavelength. 

  print Header+"Available data: "
  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = CountTheNumberOfLaserParameters(PlottingSinglePulses)
  #if(len(CEP1_uniq) == 1): #{{{ #single color simulations
  # Complexity level: 2
  
  pplt.figure(figsize=(6,10)) #new figure
  plotA=211; plotB=212
  
  ax1 = pplt.subplot(plotA)
  #ax1.set_prop_cycle('color',Colors)
  
  # Create inset in ax1 region. 
  ax11 = inset_axes(ax1, width="40%", height="30%", loc="lower right", borderpad=2)
  #ax11.set_xmargin(0.5)
  #ax11.set_ymargin(2)
  # Turning on labels for the inset
  ax11.tick_params(top=True, left=True, bottom=False, labelleft=True, labeltop=True, labelbottom=False)
  #ax11.tick_params(axis='both', top=True, left=True)
  #wpi_peaks_wavelengths  =np.zeros((len(FieldOfInterest),20)); 
  #wpi_valleys_wavelengths=np.zeros((len(FieldOfInterest),20)); 
  #wpi_peaks_Efields      =np.zeros((len(FieldOfInterest),20)); 
  #wpi_valleys_Efields    =np.zeros((len(FieldOfInterest),20)); 
  #wpi_peaks_values       =np.zeros((len(FieldOfInterest),20)); 
  #wpi_valleys_values     =np.zeros((len(FieldOfInterest),20));
  
  Peaks=[]; Valleys=[]
  Keldysh_Peaks=[]; Keldysh_Valleys=[]
  
  plot1 = []; plot11=[] #np.empty((0,0))
  for Efield1_value in Efield1_uniq: #for each value of Efield1, we make a new search, and plot one new line in the chart #{{{
    #PlottingSinglePulses_Efield = FilterDatabase(PlottingSinglePulses, Efield1_value, IndexTable['Efield1']) #select the local Efield
    PlottingSinglePulses_Efield = FilterDatabaseLowerThan(PlottingSinglePulses, Efield1_value*(1.+FieldTolerance), IndexTable['Efield1'])
    PlottingSinglePulses_Efield = FilterDatabaseGreaterThan(PlottingSinglePulses_Efield, Efield1_value*(1.-FieldTolerance), IndexTable['Efield1'])
    #print Header+"Dataset for the field: "+str(Efield1_value*1e-9)+' V/nm'
    #print PlottingSinglePulses_Efield
    # TODO: as fields are all slighly different, we should project on the neareast neighboring value, up to 1E8 V/m. 
    
    ## prints the current state of selection and stores the values of sweeped parameters
    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_Efield)
    
    # Replaces wavelengths with photon energy
    if(EnergyAxis):  #very dirty
       wavelengths1_Show = np.divide(2*np.pi*c*hbar/e, wavelengths1) #eV #E=hbar omega = 2 pi c hbar / wavelength
    else:
       wavelengths1_Show = wavelengths1
    
    ## Sorting data before plotting. 
    # If using Energies, then we sort by energy. 
    # If using wavelengths, sorting is made by wavelengths. 
    lists = sorted(itertools.izip(*[wavelengths1_Show, Ne_final, pulses1, ElectronEnergy_max]))
    wavelengths1_sorted, Ne_final_sorted, pulses1_sorted, Energy_max_sorted = list(itertools.izip(*lists))
    
    if(not EnergyAxis):
        wavelengths1_sorted_show = np.multiply(1e9, wavelengths1_sorted)
    else: #if our wavelength is actually an energy
        wavelengths1_sorted_show = wavelengths1_sorted
    wpi_sorted = np.divide(Ne_final_sorted, np.multiply(1e0, pulses1_sorted))
    
    ## Returns position and value at the peaks and valley. 
    def FindMinAndMax(signal, width): #{{{
        Header="plotFinalQuantities: FindMinAndMax():"
        print Header
        #print np.shape(signal), np.shape(width)
        width=np.asfarray(width)
        #print width
        if(EnablePeakFinder): 
          peaks, prop_peaks = find_peaks(signal) #, width)
          valleys, prop_valleys= find_peaks(-signal)
        else: 
          peaks=[]; prop_peaks=[]
          valleys=[]; prop_valleys=[]
        #print peaks, signal[peaks], width[peaks]
 
        #print "Min: "+str(valleys)
        #print "Max: "+str(peaks)
        return width[peaks], signal[peaks], width[valleys], signal[valleys]
    #}}}   
    #array_id=0
    for Efield1_compare in FieldOfInterest: #{{{ #WARNING: Efield1_value should be considered randomized
        if(abs(Efield1_value - Efield1_compare)/Efield1_value < FieldTolerance): #{{{
            wavelengths_peaks, wpi_peaks, wavelengths_valleys, wpi_valleys=FindMinAndMax(wpi_sorted, wavelengths1_sorted)
            print "W_PI peaks are at: "+str(wavelengths_peaks)
            print "W_PI hills are at: "+str(wavelengths_valleys)
            NumberOfPeaks   = np.shape(wavelengths_peaks)
            NumberOfValleys = np.shape(wavelengths_valleys)
            print "Adding extrema to a list for plotting..."
            for element in wavelengths_peaks:
                Peaks.append(  [Efield1_compare, element]) #one line per datapoint
            for element in wavelengths_valleys:
                Valleys.append([Efield1_compare, element]) #one line per datapoint
            #wpi_peaks_wavelengths  [array_id]=wavelengths_peaks  
            #wpi_peaks_values       [array_id]=wpi_peaks          
            #wpi_peaks_Efields      [array_id]=Efield1_compare    
            #wpi_valleys_wavelengths[array_id]=wavelengths_valleys
            #wpi_valleys_values     [array_id]=wpi_valleys        
            #wpi_valleys_Efields    [array_id]=Efield1_compare    
            
            #array_id+=1 #we need something to localize the values
        #### PLOTTING
            #n0 = 4.*Density
            #if(not EnergyAxis):
            plot_field, = ax1.plot(wavelengths1_sorted_show, wpi_sorted,'-^', antialiased=antialias, color=Colors[Efield1_compare], linewidth=LineWidthCoeff_abs*LineWidths[Efield1_compare], label=r'$\widetilde{w}_{\mathrm{PI}}$, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
            
            #plot_field, = ax1.plot(wavelengths1_sorted_show, wpi_sorted,'-^', linewidth=LineWidthCoeff_abs*LineWidths[Efield1_compare], label=r'$\widetilde{w}_{\mathrm{PI}}$, I='+str(FieldToIntensity(Efield1_compare))+r' W/m$^2$')
            plot1 += [plot_field]
            
            ## Let's plot here the Stark effect lines
            
            
            #plot_field, = ax1.plot(wavelengths1_Keldysh_Show, np.multiply(
                    ##PeakWpi_Gruzdev, np.divide(np.add(n0, np.multiply(-1,PeakWpi_Gruzdev*pulses1_sorted[0])), n0)),
                        #'--.'+Colors[Efield1_compare], label=r'$w_{\mathrm{PI}}\frac{n_0-n_{\mathrm{exc}}{n_0}$, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
            #plot1 += [plot_field]
            
            
            
            #plot_field, = ax1.plot(wavelengths1_Keldysh_Show, PeakWpi_Gulley, '-.'+Colors[Efield1_compare], label=r'$\widetilde{w}_{\mathrm{PI}}$ [Gulley], I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
            #plot1 += [plot_field]
            
            #plot_field, = ax1.plot(wavelengths1_Keldysh_Show, PeakWpi_Tunneling_Direct, '.-'+Colors[Efield1_compare], label=r'$\widetilde{w}_{\mathrm{TI}}^{\mathrm{(d)}}$, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
            #plot1 += [plot_field]
            #plot_field, = ax1.plot(wavelengths1_Keldysh_Show, PeakWpi_Tunneling_Indirect, '-.'+Colors[Efield1_compare], label=r'$\widetilde{w}_{\mathrm{TI}}^{\mathrm{(ind)}}$, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
            #plot1 += [plot_field]
            #plot_field, = ax1.plot(wavelengths1_Keldysh_Show, PeakWpi_Effective, '-.'+Colors[Efield1_compare], label=r'$\widetilde{w}_{\mathrm{eff}}$, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
            #plot1 += [plot_field]
            
            ## SENSELESS. 
            #plot_field, = ax1.plot(wavelengths1_Keldysh_Show, PeakWpi_Gruzdev_Indirect, '--.'+Colors[Efield1_compare], label=r'$\widetilde{w}_{\mathrm{PI}}$ [Gruzdev], I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
            #plot1 += [plot_field]
            
            #}}} #end if
            
    #}}} #end for
    del PlottingSinglePulses_Efield
  #}}} #end for
        
  Peaks=np.array(Peaks).T; Valleys=np.array(Valleys).T
  print np.shape(Peaks), np.shape(Valleys)
  ##Peaks_wavelengths=Peaks[0]
  if(PlotValleysInset):
    ax11.scatter(np.multiply(1E-9,Valleys[0,:]), Valleys[1,:], s=10, c="red", marker='^', label="Valleys TDDFT")
  else: 
    ax11.scatter(np.multiply(1E-9,Peaks[0,:]), Peaks[1,:], s=10, c="red", marker='^', label="Peaks TDDFT")
  #plot11+=[plot_field]
  #ax11.plot(wpi_peaks_Efields,   wpi_peaks_wavelengths  )
  
  #exit()
  # We separate TDDFT and Keldysh predictions now. 
  ax2 = pplt.subplot(plotB)
  # Create inset in ax2 region. 
  ax22 = inset_axes(ax2, width="30%", height="40%", loc="lower right", borderpad=2)
  
  # Turning on labels for the inset
  ax22.tick_params(labelleft=True, labelbottom=True)
  
  ## Preparing the data for Keldysh computations
  order = 50  #Keldysh series order
  
  if(PlotStarkShift): #{{{
    # Can we plot Wpi_ResonantGap(EgapEff) ? 
    IntegerLevels=np.arange(1,6,1) #generates 1,2,3,...,10
    for IntegerLevel in IntegerLevels:
        EgapResonant=Egap_SI/IntegerLevel #these are my levels resonant with initial gap
        EgapResonant_Wavelength=h*c/EgapResonant #gives resonant wavelength with field
        
        epsilon_Si0_Stark = ImportPalikDatabase_epsilon_twofiles_eV(EgapResonant_Wavelength, "Database/Si-TDLDA")
        OpticalIndex_Stark = np.real(np.sqrt(epsilon_Si0_Stark)) #maybe problem on the sqrt(real(sqrt(eps)))?
        
        # We compute Wpi_Gruzdev at this multiply resonant wavelength
        Wpi_ResonantGap = GenerateKeldyshDatabase(Egap_SI, meff_Waroquiers_LDA, EgapResonant_Wavelength, np.multiply(FieldOfInterest, np.sqrt(OpticalIndex_Stark)), order, 1.) #NOTE that OpticalIndex*FieldOfInterest was ignored for now. 
        gammaResonant = gammaKeldysh(Egap_SI, meff_Waroquiers_LDA, np.multiply(FieldOfInterest, np.sqrt(OpticalIndex_Stark)), EgapResonant_Wavelength, 1.) #valid for scalar data
        k1 = Keldysh1phi(gammaResonant); k2 = Keldysh2theta(gammaResonant) #valid
        EgapEff = EffectiveGap(Egap_SI, k1, k2)
        del k1, k2

        ## We re-inject the new resonant wavelengths into Keldysh, to compute the "real Wpi"
        #Wpi_ResonantGap = GenerateKeldyshDatabase(Egap_SI, meff_Waroquiers_LDA, h*c/EgapEff, np.multiply(FieldOfInterest, np.sqrt(OpticalIndex_Stark)), order, 1.) #NOTE that
        
        plot_field, = ax1.plot((EgapEff/e/IntegerLevel), Wpi_ResonantGap, '--k', antialiased=antialias, linewidth=1)
        plot1 += [plot_field]
        plot_field, = ax2.plot((EgapEff/e/IntegerLevel), Wpi_ResonantGap, '--k', antialiased=antialias, linewidth=1)
        plot1 += [plot_field]
        del EgapEff, Wpi_ResonantGap
  #}}}
  

  #wavelengths1_min = 200e-9; wavelengths1_max = 3300e-9 #2100e-9
  wavelengths1_Keldysh = np.linspace(wavelengths1_min, wavelengths1_max, 500)
  #wavelengths1_min = 1e-6; wavelengths1_max = 20e-6 #3300e-9 #2100e-9

  ### Capturing the dielectric permittivity in experimental or TDLDA database
  #epsilon_Si0 = ImportPalikDatabase_epsilon_fromNK(wavelengths1_Keldysh) ## TODO: replace by LDA !
  epsilon_Si0 = ImportPalikDatabase_epsilon_twofiles_eV(wavelengths1_Keldysh, "Database/Si-TDLDA")
  OpticalIndex = np.real(np.sqrt(epsilon_Si0))
    
  print Header+"Building comparison with Keldysh-Gruzdev..."
  #wavelengths1_min = np.min(wavelengths1); wavelengths1_max = np.max(wavelengths1)
  print "Dielectric permittivity of Si at 800 nm"
  print "[Palik] | [Reining et al PRL 2002]"
  print str(ImportPalikDatabase_epsilon_fromNK(800e-9))+" | "+str(ImportPalikDatabase_epsilon_twofiles_eV(800e-9, "Database/Si-TDLDA"))
    
  
  if(not AccountForOpticalIndex): 
         OpticalIndex = np.divide(OpticalIndex, OpticalIndex)
         
  for Efield1_compare in FieldOfInterest: #{{{ #INFO: good for Fig3b. 
    #NOTE: actually we can build comparison using Keldysh for much more data points. Only TDDFT limits us in number of points. And there is nothing about TDDFT in this loop! 
    
    ### I will need intermerdiate Keldysh coefficients
    #gamma = gammaKeldysh(Egap, meff, PeakField, wavelength, RefractiveIndex.real) #valid for scalar data
    #k1 = Keldysh1phi(gamma); k2 = Keldysh2theta(gamma)
    #EgapEff=EffectiveGap(Egap, k1, k2)
    
    ## We construct Wpi(corresponding to the photon energy equal to the effective gap), and plot that! 
    
    #PeakWpi_Gruzdev_Stark = 
    
    PeakWpi_Gruzdev = GenerateKeldyshDatabase(Egap_SI, EffectiveMass[Efield1_compare], wavelengths1_Keldysh, Efield1_compare*np.sqrt(OpticalIndex), order, 1., False) #NOTE: Refractive index is accounted in field, not in adiabadicity coefficient
    PeakWpi_Gruzdev_frozen = GenerateKeldyshDatabase(Egap_SI, EffectiveMass[Efield1_compare], wavelengths1_Keldysh, Efield1_compare*np.sqrt(OpticalIndex), order, 1., True) #NOTE: Refractive index is accounted in field, not in adiabadicity coefficient
    Saturation = (n0 - (PeakWpi_Gruzdev * pulseDuration )) / n0 #saturation term
    PeakWpi_Gruzdev_sat = PeakWpi_Gruzdev * Saturation
    #PeakWpi_Gulley = GenerateKeldyshGulleyDatabase(Egap_SI, EffectiveMass[Efield1_compare], wavelengths1_Keldysh, Efield1_value, order, 1.)
    #PeakWpi_Gruzdev_Indirect = GenerateKeldyshDatabase(Egap_SI_indirect, meff, wavelengths1_Keldysh, Efield1_value, order)
    
    PeakWpi_Tunneling_Direct = np.clip(KeldyshTunnelingLimit(Egap_SI, meff_Tunnel, wavelengths1_Keldysh, Efield1_compare, False, 1.), 1E36, 1E46) #NOTE: physically, refractive index does not affect tunneling time
    
    #PeakWpi_Tunneling_Indirect = np.clip(KeldyshTunnelingLimit(Egap_SI_indirect, meff_Tunnel, wavelengths1_Keldysh, Efield1_value), 1E36, 1E46)
    #PeakWpi_Gruzdev=
    #PeakWpi_Effective = np.power(np.power(PeakWpi_Gruzdev, -1) + np.power(PeakWpi_Tunneling_Direct,-1),-1)
    # Dirty: choosing energy or wavelengths. Sorting the results with this choice. 
    if(EnergyAxis):
        wavelengths1_Keldysh_Show = np.divide(2.*np.pi*c*hbar/e, wavelengths1_Keldysh) #eV
        wavelengths1_Keldysh_frozen_Show = np.divide(2.*np.pi*c*hbar/e, wavelengths1_Keldysh) #eV
    else:
        wavelengths1_Keldysh_Show = np.multiply(1E9, wavelengths1_Keldysh)
        wavelengths1_Keldysh_frozen_Show = np.multiply(1E9, wavelengths1_Keldysh)
    
    # resorting the results to plot as function of energy if necessary
    lists = sorted(itertools.izip(*[wavelengths1_Keldysh_Show, PeakWpi_Gruzdev]))
    wavelengths1_Keldysh_Show, PeakWpi_Gruzdev_Show = list(itertools.izip(*lists)) 
    
    lists_frozen = sorted(itertools.izip(*[wavelengths1_Keldysh_frozen_Show, PeakWpi_Gruzdev_frozen]))
    wavelengths1_Keldysh_frozen_Show, PeakWpi_Gruzdev_frozen_Show = list(itertools.izip(*lists_frozen)) 
    
    plot_field, = ax2.plot(wavelengths1_Keldysh_Show, PeakWpi_Gruzdev_Show, '-', antialiased=antialias, color=Colors[Efield1_compare], linewidth=LineWidthCoeff_rel*LineWidthCoeff_abs*LineWidths[Efield1_compare], label=r'$\widetilde{w}_{\mathrm{PI}}^{\mathrm{KG}}$, I='+str(FieldToIntensity(Efield1_compare))+r' W/cm$^2$')
    plot1 += [plot_field]
    
    ## Plotting Keldysh-Gruzdev theory in inset of Fig (bottom) (ax22). 
    plot_field, = ax22.plot(wavelengths1_Keldysh_frozen_Show, PeakWpi_Gruzdev_frozen_Show, '-', antialiased=antialias, color=Colors[Efield1_compare], linewidth=LineWidthCoeff_rel*LineWidthCoeff_abs*LineWidths[Efield1_compare], label=r'$\widetilde{w}_{\mathrm{PI}}^{\mathrm{KG}}$ frozen, I='+str(FieldToIntensity(Efield1_compare))+r' W/cm$^2$')

    #plot_field, = ax1.plot(wavelengths1_Keldysh_Show, PeakWpi_Gruzdev_Show, '-', linewidth=LineWidthCoeff_rel*LineWidthCoeff_abs*LineWidths[Efield1_compare], label=r'$\widetilde{w}_{\mathrm{PI}}$ [Gruzdev], I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
    
    plot1 += [plot_field]
  #}}} #end for plotting Fig2b
  
  #INFO: provides the Keldysh analysis part for Fig3a (inset).
  #WARNING: here we have disabled the contribution of optical index. Then one obtain the match for gap resonances with TDDFT. 
  # However, it seems W_PI amplitude is good only if we account for optical index. How the optical index should be treated in Keldysh theory or in TDDFT ?
  for Efield1_compare in np.linspace(0,np.max(FieldOfInterest),100):
    
    PeakWpi_Gruzdev_continuous = GenerateKeldyshDatabase(Egap_SI, meff_Waroquiers_LDA, wavelengths1_Keldysh, Efield1_compare*np.sqrt(OpticalIndex), order, 1., False) #NOTE: Refractive index is accounted in field, not in adiabadicity coefficient
    
    # Dirty: choosing energy or wavelengths. Sorting the results with this choice. 
    if(EnergyAxis):
        wavelengths1_Keldysh_Show = np.divide(2.*np.pi*c*hbar/e, wavelengths1_Keldysh) #eV
        wavelengths1_Keldysh_frozen_Show = np.divide(2.*np.pi*c*hbar/e, wavelengths1_Keldysh) #eV
    else:
        wavelengths1_Keldysh_Show = np.multiply(1E9, wavelengths1_Keldysh)
        wavelengths1_Keldysh_frozen_Show = np.multiply(1E9, wavelengths1_Keldysh)
    
    # resorting the results to plot as function of energy if necessary
    lists = sorted(itertools.izip(*[wavelengths1_Keldysh_Show, PeakWpi_Gruzdev_continuous]))
    wavelengths1_Keldysh_Show, PeakWpi_Gruzdev_Show = list(itertools.izip(*lists)) 
    
    lists_frozen = sorted(itertools.izip(*[wavelengths1_Keldysh_frozen_Show, PeakWpi_Gruzdev_frozen]))
    wavelengths1_Keldysh_frozen_Show, PeakWpi_Gruzdev_frozen_Show = list(itertools.izip(*lists_frozen)) 
    
    print np.shape(PeakWpi_Gruzdev_Show), np.shape(wavelengths1_Keldysh_Show)
    wavelengths_peaks, wpi_peaks, wavelengths_valleys, wpi_valleys=FindMinAndMax(np.array(PeakWpi_Gruzdev_Show), np.array(wavelengths1_Keldysh_Show))
    
    print "W_PI peaks are at: "+str(wavelengths_peaks)
    print "W_PI hills are at: "+str(wavelengths_valleys)
    NumberOfPeaks   = np.shape(wavelengths_peaks)
    NumberOfValleys = np.shape(wavelengths_valleys)
    print "Adding extrema to a list for plotting..."
    for element in wavelengths_peaks:
        Keldysh_Peaks.append(  [Efield1_compare, element]) #one line per datapoint
    for element in wavelengths_valleys:
        Keldysh_Valleys.append([Efield1_compare, element]) #one line per datapoint
    
    
    plot1 += [plot_field]
  #}}} #end for plotting Fig1a inset.
  
  Keldysh_Peaks=np.array(Keldysh_Peaks).T; Keldysh_Valleys=np.array(Keldysh_Valleys).T
  print np.shape(Keldysh_Peaks), np.shape(Keldysh_Valleys)
  
  if(PlotValleysInset):
      ax11.scatter(np.multiply(1E-9,Keldysh_Valleys[0,:]), Keldysh_Valleys[1,:], s=0.3, c="blue", marker='*', label="Valleys Keldysh")
  else:
      ax11.scatter(np.multiply(1E-9,Keldysh_Peaks[0,:]), Keldysh_Peaks[1,:], s=0.3, c="blue", marker='*', label="Peaks Keldysh")
    
  #if(EnergyAxis): 
    #ax1.set_xscale('log')
  ax1.set_yscale('log')
  ax1.set_ylabel(r'$\tilde{w}_{PI}$ (m$^{-3}$.s$^{-1}$)')
  ax11.set_xlabel(r'$E_{\mathrm{field}}$ (V/nm)')
  ax11.set_ylabel(r'$\hbar \omega$ (eV)')
  #if(PlotValleysInset):
      #ax11.set_title(r"Min of $\tilde{w}_{PI}$")
  #else:
      #ax11.set_title(r"Max of $\tilde{w}_{PI}$")
  ##ax11.set_xscale('log')
  #ax11.set_yscale('log')
  ax2.set_ylabel(r'$\tilde{w}_{PI}$ (m$^{-3}$.s$^{-1}$)')
  ax2.set_yscale('log')
  ax22.set_yscale('log')
  
  if(ShowDensities): #{{{
    ax12 = ax1.twinx()
    #pplt.title('Final density as function of laser wavelength')
    print Header+"Total number of datasets: "+str(len(set(Efield1_uniq)))
    
    for Efield1_value in Efield1_uniq: #for each value of Efield1, we make a new search, and plot one new line in the chart
        #PlottingSinglePulses_Efield = FilterDatabase(PlottingSinglePulses, Efield1_value, IndexTable['Efield1']) #select the local Efield
        # TODO: how to merge data with values close to each others? Maybe introduce tolerance of 10% ?
        PlottingSinglePulses_Efield = FilterDatabaseLowerThan(PlottingSinglePulses, Efield1_value*(1.+FieldTolerance), IndexTable['Efield1'])
        PlottingSinglePulses_Efield = FilterDatabaseGreaterThan(PlottingSinglePulses_Efield, Efield1_value*(1.-FieldTolerance), IndexTable['Efield1'])
        
        #print PlottingSinglePulses_Efield
        ## prints the current state of selection and stores the values of sweeped parameters
        runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_Efield)
        
        if(EnergyAxis): #very dirty!
          wavelengths1=np.divide(2.*np.pi*c*hbar/e, wavelengths1)
        
        lists = sorted(itertools.izip(*[wavelengths1, Ne_final, ElectronEnergy_max])) #NOTE: Ne_max makes more sense than Ne_final
        wavelengths1_sorted, Ne_final_sorted, Energy_max_sorted = list(itertools.izip(*lists))
        if(not EnergyAxis): 
            wavelengths1_sorted = np.multiply(1e9, wavelengths1_sorted)
        
        ax12.set_ylabel(r'Final density $n_{\mathrm{exc}}$ (m$^{-3}$)')
        #print Efield1_value
        # TODO: PUT here a select case which attributes a new color to a new field interval
        # OR USE A DICTIONNARY. 
        
        # SAME CODE, just more clever. 
        for Efield1_compare in FieldOfInterest:
            if(abs(Efield1_value - Efield1_compare)/Efield1_value < FieldTolerance):
                #if(not EnergyAxis):
                plot_field, = ax12.plot(wavelengths1_sorted, Ne_final_sorted, '^', color=Colors[Efield1_compare], label=r'$E_{peak}=$'+str(round(Efield1_compare*1E-9))+' V/nm')

                #plot_field, = ax12.plot(wavelengths1_sorted, Ne_final_sorted, '^', label=r'$E_{peak}=$'+str(round(Efield1_compare*1E-9))+' V/nm')
                plot1 += [plot_field]
                #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
                break
        
        del PlottingSinglePulses_Efield
  #}}}  
  labels1 = [l.get_label() for l in plot1]
  #ax1.legend(plot1, labels1, loc='best') #outside right, middle
  #ax1.legend(loc=3, borderaxespad=0., bbox_to_anchor=(1.05, 0.05))
  if(not EnergyAxis): 
    ax1.set_xlabel('Laser wavelength (nm)')
    ax1.set_xlim([wavelengths1_min, wavelengths1_max])
    ax11.set_ylim([wavelengths1_min,wavelengths1_max])
    ax2.set_xlabel('Laser wavelength (nm)')
    ax2.set_xlim([wavelengths1_min, wavelengths1_max])
    ax22.set_xlim([wavelengths1_min,wavelengths1_max])
  else: 
    ax1.set_xlabel('Photon energy (eV)')
    ax1.set_xlim([h*c/wavelengths1_max/e, h*c/wavelengths1_min/e])
    #ax11.set_ylim([h*c/wavelengths1_max/e, h*c/wavelengths1_min/e])
    ax11.set_ylim([h*c/wavelengths1_max/e, 1.5])
    ax2.set_xlabel('Photon energy (eV)')
    ax2.set_xlim([h*c/wavelengths1_max/e, h*c/wavelengths1_min/e])
    ax22.set_xlim([h*c/wavelengths1_max/e, h*c/wavelengths1_min/e])
  ax1.set_ylim([1E35,1E43])
  ax2.set_ylim([1E35,1E43])
  ax22.set_ylim([1E35,1E43])
  #ax1.set_xscale('log')
  #ax1.set_yscale('log')
  #ax1.grid()
  
  ResonantWavelength = h*c/Egap_SI
  print "** Adding vertical bar for the resonant transitions at: "+str(round(1e9*ResonantWavelength))+"nm."
  
  if(not EnergyAxis):
    ax1.axvline(x=1e9*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax1.axvline(x=1e9*2*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax1.axvline(x=1e9*3*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax1.axvline(x=1e9*4*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    
    ax2.axvline(x=1e9*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax2.axvline(x=1e9*2*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax2.axvline(x=1e9*3*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax2.axvline(x=1e9*4*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    
    ax22.axvline(x=1e9*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax22.axvline(x=1e9*2*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax22.axvline(x=1e9*3*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax22.axvline(x=1e9*4*ResonantWavelength, color='k', linestyle='--', linewidth=1)
  else:
    ax1.axvline(x=Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax1.axvline(x=1./2.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax1.axvline(x=1./3.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax1.axvline(x=1./4.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax11.axhline(y=Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax11.axhline(y=1./2.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax11.axhline(y=1./3.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax11.axhline(y=1./4.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    
    ax2.axvline(x=Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax2.axvline(x=1./2.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax2.axvline(x=1./3.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax2.axvline(x=1./4.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    
    ax22.axvline(x=Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax22.axvline(x=1./2.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax22.axvline(x=1./3.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax22.axvline(x=1./4.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
  
  ## Saving this first double figure
  pplt.tight_layout()
  filename = 'Metaplot-SinglePulse-SingleColor-EffectOfWavelength_Tol10-Keldysh-'+str(int(pulseDuration*1E15))
  if(PlotValleysInset):
      filename+="-Valleys"
  else:
      filename+="-Hills"
  if(AccountForOpticalIndex): 
      filename+="-WithOpticalIndex"
  else:
      filename+="-NoOpticalIndex"
  pplt.savefig(filename+'.eps')
  pplt.savefig(filename+'.png')
  
  ## Imports the optical data from Palik database
  # Returns the spectrum (in Re(eps), Im(eps)) of a material on the given grid of wavelengths

  ## Computes the free carrier energy contained in the volume <SampleVolume>.  
  def EnergyFromFreeCarrierAbsorption(pulseDuration, wavelengths1_sorted, Efield, Ne_final_sorted, w_PI_sorted, Egap_SI = 2.56*e, meff = 0.2226, SampleVolume=2 * 5E28**-1): #{{{
    # dEnergy/dt = alpha * intensity + ( n har omega - Egap ) x sigma_n I^n
    # where we approach sigma_n I^n by pulse-averaged w_PI 
    # Energy ~ alpha * intensity + ( n har omega - Egap ) x w_PI * tau
    omega   = np.divide(2. * np.pi * c, wavelengths1_sorted)
    #print Header+"omega: "+str(omega) #Valid
    nu = 6e-15**-1 #1.1 fs: Solowski-Tinten adjustment. Can also use the one from Shunsuke Sato, PRB 2015 (Numerical pump probe)
    print "nu: "+str(nu) #Valid
    #epsilon_Si0 = 13.64 + 0.048j #TODO: get the optical data in extended spp theory for Si
    epsilon_Si0 = ImportPalikDatabase_epsilon_fromNK(wavelengths1_sorted) ## TODO: replace by LDA !
    print "Dielectric permittivity of Si at 800 nm : "+str(ImportPalikDatabase_epsilon_fromNK(800e-9))
    omega_p2    = np.multiply(Ne_final_sorted, e**2 / m_e / meff / epsilon_0)
    #print "omega_p: "+str(np.sqrt(omega_p2)) # Valid
    # NOTE: THIS IS A CHOICE. WE JUST HAVE TO BE CONSISTENT ON WHICH PROCESS WE WISH FOR. 
    Epsilon = 1E0 - np.divide(omega_p2, np.power(omega, 2)) * np.divide(1.,1. + np.divide(1.j * nu, omega))
    #Epsilon = epsilon_Si0 - np.divide(omega_p2, np.power(omega, 2)) * np.divide(1.,1. + np.divide(1.j * nu, omega))
    #NumberOfPhotons = Egap_SI / (hbar * omega)
    NumberOfPhotons = np.divide(Egap_SI/hbar, omega)
    print NumberOfPhotons
    #NumberOfPhotons_round = np.round(NumberOfPhotons,0)
    NumberOfPhotons_round = np.ceil(NumberOfPhotons)
    print Header+"sqrt(Epsilon): "+str(np.min(np.sqrt(Epsilon)))+", "+str(np.max(np.sqrt(Epsilon))) #Valid
    alpha                = np.multiply(2.* omega/c, np.imag(np.sqrt(Epsilon))) #valid
    EnergyExcess         = np.multiply(NumberOfPhotons_round * hbar * omega, Ne_final_sorted) * SampleVolume / e  #result in eV
    EnergyBremsstrahlung = np.multiply( omega * epsilon_0 * (Epsilon.imag), np.multiply(Efield, Efield)) * SampleVolume / e  * pulseDuration
    EnergyBremsstrahlung_approx =  EnergyBremsstrahlung #NOTE: not used for now. Another collision freq? 
    #EnergyBremsstrahlung_approx = np.multiply( 4.*np.pi / c * np.divide(Epsilon.imag, (np.sqrt(Epsilon)).real), )) ) * SampleVolume / e
    print Header+"FcrEnergy (eV): "+str(np.min(EnergyBremsstrahlung))+", "+str(np.max(EnergyBremsstrahlung))
    print Header+"Excess energy (n hbar omega) (eV): "+str(np.min(EnergyExcess))+", "+str(np.max(EnergyExcess))
    EnergyFreeCarrier    = EnergyExcess + EnergyBremsstrahlung
    EnergyFreeCarrier_approx = EnergyExcess + EnergyBremsstrahlung_approx
    return EnergyFreeCarrier, EnergyFreeCarrier_approx
  #}}}
  
  
  
  
  ## TODO: TDDFT vs Drude should become another figure. We don't need to have them side by side. 
  plt.figure()
  ax3 = plt.subplot(111)
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
    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_Efield)
    
    # 
    if(EnergyAxis):
        wavelengths1_Show = np.divide(2.*np.pi*c*hbar/e, wavelengths1) #eV
    else: 
        wavelengths1_Show = np.multiply(1e9, wavelengths1)
    
    lists = sorted(itertools.izip(*[wavelengths1_Show, wavelengths1, Ne_final, W_pi, ElectronEnergy_max])) #NOTE: Ne_max makes more sense than Ne_final
    wavelengths1_sorted_show, wavelengths1_sorted, Ne_final_sorted, W_pi_sorted, Energy_max_sorted = list(itertools.izip(*lists))
    #print Efield1_value
    
    for Efield1_compare in FieldOfInterest: 
        if(abs(Efield1_value - Efield1_compare)/Efield1_value < FieldTolerance):
            
            ax3.plot(wavelengths1_sorted_show, Energy_max_sorted, '-', linewidth=LineWidthCoeff_abs*LineWidths[Efield1_compare], color=Colors[Efield1_compare])
            #ax3.plot(wavelengths1_sorted_show, Energy_max_sorted, '-')
            
    del PlottingSinglePulses_Efield
    
  #pplt.xlabel(r'Intensity (W/m$^2$)')
  if(not EnergyAxis):
    ax3.set_xlabel('Laser wavelength (nm)')
  else: 
    ax3.set_xlabel('Photon energy (eV)')
  #pplt.ylabel('Irreversibility factor')
  ax3.set_ylabel('Electron excess energy (eV)')
  #pplt.ylim((0,1))
  #ax2.set_xscale('log')
  
  ## TODO: this can be removed, since irreversibility is gauge-dependent. 
  #if(PlotIrreversibility): #{{{
    #ax22 = ax2.twinx()
    #for Efield1_value in Efield1_uniq: #for each value of Efield1, we make a new search, and plot one new line in the chart
        ##PlottingSinglePulses_Efield = FilterDatabase(PlottingSinglePulses, Efield1_value, IndexTable['Efield1']) #select the local Efield
        #PlottingSinglePulses_Efield = FilterDatabaseLowerThan(PlottingSinglePulses, Efield1_value*(1.+FieldTolerance), IndexTable['Efield1'])
        #PlottingSinglePulses_Efield = FilterDatabaseGreaterThan(PlottingSinglePulses_Efield, Efield1_value*(1.-FieldTolerance), IndexTable['Efield1'])
        ##print Header+"Dataset for the field: "+str(Efield1_value*1e-9)+' V/nm'
        ##print PlottingSinglePulses_Efield
        ## TODO: as fields are all slighly different, we should project on the neareast neighboring value, up to 1E8 V/m. 
        
        ### prints the current state of selection and stores the values of sweeped parameters
        #runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, ElectronEnergy_max, Energy_max, W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_Efield)
        
        #lists = sorted(itertools.izip(*[wavelengths1, Ne_final, Energy_max, irreversibility])) #NOTE: Ne_max makes more sense than Ne_final
        #wavelengths1_sorted, Ne_final_sorted, Energy_max_sorted, irreversibility_sorted = list(itertools.izip(*lists))
        #wavelengths1_sorted_show = np.multiply(1e9, wavelengths1_sorted)
        
        #for Efield_compare in FieldOfInterest:
            #if(abs(Efield1_value - Efield1_compare)/Efield1_compare < FieldTolerance):
                #ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+'+Colors[Efield1_compare], label='Irreversibility')
        
        ###print Efield1_value
        ### TODO: PUT here a select case which attributes a new color to a new field interval
        ###if(1e8*(1.-FieldTolerance) < Efield1_value < 1e8*(1.+FieldTolerance)):
            ###ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+r', label='Irreversibility')
            ###ax22.plot(wavelengths1_sorted_show, gammaKeldysh(Egap_SI, 0.228, Efield1_value, wavelengths1_sorted), '--r', label='Adiabaticity')
        ##if(1e9*(1.-FieldTolerance) < Efield1_value < 1e9*(1.+FieldTolerance)):
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+b', label='Irreversibility')
            ###ax22.plot(wavelengths1_sorted_show, gammaKeldysh(Egap_SI, 0.228, Efield1_value, wavelengths1_sorted), '--b', label='Adiabaticity')
        ##if(5e9*(1.-FieldTolerance) < Efield1_value < 5e9*(1.+FieldTolerance)):
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+g', label='Irreversibility')
            ###ax22.plot(wavelengths1_sorted_show, gammaKeldysh(Egap_SI, 0.228, Efield1_value, wavelengths1_sorted), '--g', label='Adiabaticity')
        ##if(1e10*(1.-FieldTolerance) < Efield1_value < 1e10*(1.+FieldTolerance)):
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+k', label='Irreversibility')
            ###ax22.plot(wavelengths1_sorted_show, gammaKeldysh(Egap_SI, 0.228, Efield1_value, wavelengths1_sorted), '--k', label='Adiabaticity')
        #del PlottingSinglePulses_Efield
    
    
    #pplt.ylabel(r'Irreversibility $\zeta$')
    #ax22.set_yscale('log')
  ##}}}
  

  if(PlotFreeCarrierAbsorption):
      print Header+"Now plotting the free carrier energy absorption..."
      for Efield1_value in Efield1_uniq: #for each value of Efield1, we make a new search, and plot one new line in the chart
        #PlottingSinglePulses_Efield = FilterDatabase(PlottingSinglePulses, Efield1_value, IndexTable['Efield1']) #select the local Efield
        PlottingSinglePulses_Efield = FilterDatabaseLowerThan(PlottingSinglePulses, Efield1_value*(1.+FieldTolerance), IndexTable['Efield1'])
        PlottingSinglePulses_Efield = FilterDatabaseGreaterThan(PlottingSinglePulses_Efield, Efield1_value*(1.-FieldTolerance), IndexTable['Efield1'])
        #print Header+"Dataset for the field: "+str(Efield1_value*1e-9)+' V/nm'
        #print PlottingSinglePulses_Efield
        # TODO: as fields are all slighly different, we should project on the neareast neighboring value, up to 1E8 V/m. 
        
        ## prints the current state of selection and stores the values of sweeped parameters
        runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_Efield)
        
        if(EnergyAxis):
            wavelengths1_Show = np.divide(2.*np.pi*c*hbar/e, wavelengths1) #eV
        else: 
            wavelengths1_Show = np.multiply(1e9, wavelengths1)
    
        lists = sorted(itertools.izip(*[wavelengths1_Show, wavelengths1, Ne_final, W_pi, ElectronEnergy_max])) #NOTE: Ne_max makes more sense than Ne_final
        wavelengths1_sorted_show, wavelengths1_sorted, Ne_final_sorted, W_pi_sorted, Energy_max_sorted = list(itertools.izip(*lists))
    #print Efield1_value
    
        
        #lists = sorted(itertools.izip(*[wavelengths1, Ne_final, W_pi, Energy_max, irreversibility])) #NOTE: Ne_max makes more sense than Ne_final
        #wavelengths1_sorted, Ne_final_sorted, W_pi_sorted, Energy_max_sorted, irreversibility_sorted = list(itertools.izip(*lists))
        #wavelengths1_sorted_show = np.multiply(1e9, wavelengths1_sorted)
        
        for Efield1_compare in FieldOfInterest: 
            if(abs(Efield1_value - Efield1_compare) / Efield1_value < FieldTolerance):
                EnergyAnalytical, EnergyAnalytical_approx = EnergyFromFreeCarrierAbsorption(pulses1, wavelengths1_sorted, Efield1_value, Ne_final_sorted, W_pi_sorted, Egap_SI, EffectiveMass[Efield1_compare])
                                
                ax3.plot(wavelengths1_sorted_show, EnergyAnalytical, '--', color=Colors[Efield1_compare], label='Free carrier absorption')
                #ax2.plot(wavelengths1_sorted_show, EnergyAnalytical, '--', label='Free carrier absorption')
            
        ##if(1e8*(1.-FieldTolerance) < Efield1_value < 1e8*(1.+FieldTolerance)):
            ##print Header+"Found field 0.1 V/nm"
            ##EnergyAnalytical, EnergyAnalytical_approx = EnergyFromFreeCarrierAbsorption(pulses1, wavelengths1_sorted, Efield1_value, Ne_final_sorted, W_pi_sorted, Egap_SI)
            ##ax2.plot(wavelengths1_sorted_show, EnergyAnalytical, '--+r', label='Free carrier absorption')
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+r', label='Irreversibility')
            ##ax22.plot(wavelengths1_sorted_show, gammaKeldysh(Egap_SI, 0.228, Efield1_value, wavelengths1_sorted), '--r', label='Adiabaticity')
        #if(1e9*(1.-FieldTolerance) < Efield1_value < 1e9*(1.+FieldTolerance)):
            #print Header+"Found field 1.0 V/nm"
            #EnergyAnalytical, EnergyAnalytical_approx = EnergyFromFreeCarrierAbsorption(pulses1, wavelengths1_sorted, Efield1_value, Ne_final_sorted, W_pi_sorted, Egap_SI)
            #ax2.plot(wavelengths1_sorted_show, EnergyAnalytical, '--+b', label='Free carrier absorption')
            ##ax2.plot(wavelengths1_sorted_show, EnergyAnalytical_approx, '--b', label='Free carrier absorption')
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+b', label='Irreversibility')
            ##ax22.plot(wavelengths1_sorted_show, gammaKeldysh(Egap_SI, 0.228, Efield1_value, wavelengths1_sorted), '--b', label='Adiabaticity')
        #if(5e9*(1.-FieldTolerance) < Efield1_value < 5e9*(1.+FieldTolerance)):
            #print Header+"Found field 5 V/nm"
            #EnergyAnalytical, EnergyAnalytical_approx = EnergyFromFreeCarrierAbsorption(pulses1, wavelengths1_sorted, Efield1_value, Ne_final_sorted, W_pi_sorted, Egap_SI)
            #ax2.plot(wavelengths1_sorted_show, EnergyAnalytical, '--+g', label='Free carrier absorption')
            ##ax2.plot(wavelengths1_sorted_show, EnergyAnalytical_approx, '--g', label='Free carrier absorption')
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+g', label='Irreversibility')
            ##ax22.plot(wavelengths1_sorted_show, gammaKeldysh(Egap_SI, 0.228, Efield1_value, wavelengths1_sorted), '--g', label='Adiabaticity')
        #if(1e10*(1.-FieldTolerance) < Efield1_value < 1e10*(1.+FieldTolerance)):
            #print Header+"Found field 10 V/nm"
            #EnergyAnalytical, EnergyAnalytical_approx = EnergyFromFreeCarrierAbsorption(pulses1, wavelengths1_sorted, Efield1_value, Ne_final_sorted, W_pi_sorted, Egap_SI)
            #ax2.plot(wavelengths1_sorted_show, EnergyAnalytical, '--+k', label='Free carrier absorption')
            ##ax2.plot(wavelengths1_sorted_show, EnergyAnalytical_approx, '--k', label='Free carrier absorption')
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+k', label='Irreversibility')
            ##ax22.plot(wavelengths1_sorted_show, gammaKeldysh(Egap_SI, 0.228, Efield1_value, wavelengths1_sorted), '--k', label='Adiabaticity')
        del PlottingSinglePulses_Efield
    
    #pplt.ylabel(r'Irreversibility $\zeta$')
    #ax22.set_yscale('log')
      #ax22 = ax2.twinx()
  #pplt.ylabel(r'Adiabaticity $\gamma$, Irreversibility $\zeta$')
  #ax22.set_ylim((0,1))
  ax3.set_xscale("log")
  ax3.set_yscale("symlog", linthreshy=1e-9)
  #ax2.set_ylim((1E-9,5E2))
  #ax22.grid()
  if(not EnergyAxis): 
    pplt.xlim((200.,2100.))
  else: 
    pplt.xlim((0.3,10.))
  
  
  #ResonantWavelength = h*c/Egap_SI
  print "** Adding vertical bar for the resonant transitions at: "+str(round(1e9*ResonantWavelength))+"nm."
  
  if(not EnergyAxis):
    ax3.axvline(x=1e9*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax3.axvline(x=1e9*2*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax3.axvline(x=1e9*3*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax3.axvline(x=1e9*4*ResonantWavelength, color='k', linestyle='--', linewidth=1)
  else:
    ax3.axvline(x=Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax3.axvline(x=1./2.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax3.axvline(x=1./3.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax3.axvline(x=1./4.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
  #ax2.legend(loc=2, borderaxespad=0.)
  #ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
  #ax22.legend(loc=1, borderaxespad=0.)

  # Saving the plot
  pplt.tight_layout
  filename = 'Metaplot-SinglePulse-SingleColor-EffectOfWavelength_Tol10-Drude-'+str(int(pulseDuration*1E15))
  pplt.savefig(filename+'.eps')
  pplt.savefig(filename+'.png')
  if(ShowFinalPlot==1):
    pplt.show()

  del ax1, PlottingSinglePulses
  #}}}

## Plots (Wpi | Energy) (wavelength) (1D plot) for the special case of 2-color irradiations. Fields are selected by introducing <FieldToTolerance> percents tolerance on the electric field amplitude. 
# This enables to analyze the role of laser wavelength on the final excited electron density, excitation rate and laser-induced energy increase. Provides predictions on the Stark effect. 
# Energy is compared with analytical Drude model. 
def Bicolor_plotEffectOfWavelength_Tolerance(PlottingTable, wavelength1, FieldTolerance=0.05, Field1=1E0, Field2max=1E9, pulseDuration=20e-15, CEP1=0., PulseType=3, RiseTime=2, ShowFinalPlot=0, Egap_SI=2.56*e, Egap_SI_indirect=0.51*e): #{{{
  ShowDensities            =False
  PlotIrreversibility      =False
  PlotFreeCarrierAbsorption=False
  EnergyAxis               =False #True: Xaxis is in eV. False: wavelength in nm. 
  
  CutOffWavelength=2.1e-6
  
  CEP2_request=-np.pi/2. #np.pi/2.
  
  meff_Waroquiers_LDA = 0.2226 #effective mass of electron [Waroquiers, D. et al, Physical Review B, 2013, 87, 075121]
  #meff = 0.14 #0.14: fits for Si in MPI regime at 800 nm.
  #Field dependent fitting of effective mass
  meff_LowField = meff_Waroquiers_LDA
  meff_HighField = meff_Waroquiers_LDA #0.14: fits for Si in MPI regime at 800 nm.
  meff_VeryHighField = meff_Waroquiers_LDA
  
  meff_Tunnel = meff_Waroquiers_LDA #fits well if using INDIRECT band gap for 800 nm
  print "#====== Plot: DOUBLE PULSE, EFFECT OF THE WAVELENGTH, SAME CEP, SAME POLARIZATIONS: final density as function of field ======"
  #Filter: double pulses, same wavelengths, same polarization
  
  LineWidthCoeff_abs = 0.5 #controls absolute thickness of all curves
  LineWidthCoeff_rel = 1.0 #controls thickness of Keldysh-Gruzdev curves relatively to TDDFT
  antialias=True
  
  #FieldOfInterest = [5e8, 1e9, 2e9, 3E9, 4E9, 5e9, 1e10]
  FieldOfInterest = [2.5E8, 5E8, 1E9, 5E9] #[0E0, 1.0E9] #, 1e10] # NOTE: this should be used only for Efield2
  NumberOfColors=len(FieldOfInterest)
  ColorList=[plt.cm.seismic(i) for i in np.linspace(0, 1, NumberOfColors)] 
  # Not bad! RdYlBu, coolwarm, seismic
  # Too soft: Copper, Spectral
  # Not perfect. RdBu. RdYlGn.
  # Too dynamic: bwr
  # To test: , brg, rainbow
  #Colors={5e8: 'c', 7.5e8: 'r', 1e9: 'b', 1.5e9: 'g', 2e9: 'y', 3e9: 'm'} #, 1e10: 'k'}
  Colors={FieldOfInterest[i]:ColorList[i] for i in range(0,NumberOfColors)}
  LineWidths={FieldOfInterest[i]:np.arange(1.5,4,NumberOfColors) for i in range(0,NumberOfColors)}
  #EffectiveMass={1e9: meff_LowField, 2e9: meff_LowField, 5e9: meff_HighField, 1e10: meff_VeryHighField}
  EffectiveMass={FieldOfInterest[i]:meff_HighField for i in range(0, NumberOfColors)}
  # Simple pulse, single color
  
  # TODO: Filtering function
  PlottingDoublePulses = PlottingTable #get the full database
  
  # The wavelength1 parameter enables to choose which wavelength1 we want to observe the Stark effect for. Of course, we shall not discernate laser1 and laser2 in the data. But how? Lets discernate first. And invert later. 
  PlottingDoublePulses = FilterDatabase(PlottingDoublePulses, wavelength1, IndexTable['wavelengths1'])
  
  PlottingDoublePulses = FilterDatabaseGreaterThan(PlottingDoublePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. # NOTE: works
  
  # NOTE: 0 should be kept here, since we use non-strict inequality
  PlottingDoublePulses = FilterDatabaseGreaterThan(PlottingDoublePulses, Field1*(1-FieldTolerance), IndexTable['Efield1']) #select fields that are close to selected Efield1
  PlottingDoublePulses = FilterDatabaseLowerThan(PlottingDoublePulses, Field1*(1+FieldTolerance), IndexTable['Efield1']) #select fields that are close to selected Efield1
  
  # We keep Efield2==0 possible. 
  PlottingDoublePulses = FilterDatabaseLowerThan(PlottingDoublePulses, Field2max, IndexTable['Efield2']) #we want to include cases where Efield2==0, and Efield2 >> 1. However, we want to put a maximum field for Efield2. 
  
  PlottingDoublePulses = FilterDatabase(PlottingDoublePulses, 0., IndexTable['CrossPolarization']) #removes fields where CrossPolarization is different from 0. # TODO: CrossPolarization=1 case would be interesting for Stark effect. 
  
  PlottingDoublePulses = FilterDatabase(PlottingDoublePulses, PulseType, IndexTable['PulseType']) #keeps only pulses of type StH (3)
  
  PlottingDoublePulses = FilterDatabase(PlottingDoublePulses, pulseDuration, IndexTable['pulses1']) #keeps only pulses of selected duration in argument
  
  #if(RiseTime >= 1): #select by RiseTime only if it is defined
    #PlottingDoublePulses = FilterDatabase(PlottingDoublePulses, RiseTime, IndexTable['RiseTime'])
  
  PlottingDoublePulses = FilterDatabase(PlottingDoublePulses, CEP1, IndexTable['CEP1']) #keeps only CEP = 0. 
  
  if(CEP2_request >= 0):
      PlottingDoublePulses = FilterDatabaseLowerThan(PlottingDoublePulses,   (CEP2_request)*(1.1), IndexTable['CEP2'])
      PlottingDoublePulses = FilterDatabaseGreaterThan(PlottingDoublePulses, (CEP2_request)*(0.9), IndexTable['CEP2']) #keeps only CEP = 0. 
  else:
      PlottingDoublePulses = FilterDatabaseGreaterThan(PlottingDoublePulses, (CEP2_request)*(1.1), IndexTable['CEP2'])
      PlottingDoublePulses = FilterDatabaseLowerThan(PlottingDoublePulses,   (CEP2_request)*(0.9), IndexTable['CEP2']) #keeps only CEP = 0. 
  # TODO: we shall clean CEP2 as well, as it must influence the results. 
  
  PlottingDoublePulses = FilterDatabaseLowerThan(PlottingDoublePulses, CutOffWavelength, IndexTable['wavelengths2']) #introducing a cutoff at 2.1 um in selection. 

  ## TODO: check if CEP could count at smaller wavelength. 

  print Header+"Available data: "
  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = CountTheNumberOfLaserParameters(PlottingDoublePulses)
  #if(len(CEP1_uniq) == 1): #{{{ #single color simulations
  # Complexity level: 2
  
  #print("Number of Efield2 uniq: ")
  #print np.shape(Efield2_uniq)
  
  print(CEP2_uniq)
  #exit()
  
  pplt.figure(figsize=(12,6)) #new figure
  plotA=121; plotB=122
  
  ax1 = pplt.subplot(111)
  
  #ax1.set_prop_cycle('color',Colors)
    
  plot1 = [] #np.empty((0,0))
  for Efield2_compare in FieldOfInterest: #{{{
    PlottingDoublePulses_Efield = PlottingDoublePulses
    PlottingDoublePulses_Efield = FilterDatabaseLowerThan(PlottingDoublePulses_Efield, Efield2_compare*(1.+FieldTolerance), IndexTable['Efield2'])
    PlottingDoublePulses_Efield = FilterDatabaseGreaterThan(PlottingDoublePulses_Efield, Efield2_compare*(1.-FieldTolerance), IndexTable['Efield2'])
    #print Header+"Dataset for the field: "+str(Efield2_value*1e-9)+' V/nm'
    #print PlottingDoublePulses_Efield
    # TODO: as fields are all slighly different, we should project on the neareast neighboring value, up to 1E8 V/m. 

    #print Header+"Available data for the specific field E2="+str(Efield2_value*1E-9)+" V/nm."
    #tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9 = CountTheNumberOfLaserParameters(PlottingDoublePulses_Efield)

    ## prints the current state of selection and stores the values of sweeped parameters
    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingDoublePulses_Efield)
    
    # Replaces wavelengths with photon energy
    if(EnergyAxis):  #very dirty
       wavelengths2_Show = np.divide(2*np.pi*c*hbar/e, wavelengths2) #eV #E=hbar omega = 2 pi c hbar / wavelength
    else:
       wavelengths2_Show = wavelengths2
    
    ## Sorting data before plotting. 
    # If using Energies, then we sort by energy. 
    # If using wavelengths, sorting is made by wavelengths. 
    lists = sorted(itertools.izip(*[wavelengths2_Show, Ne_final, pulses1, ElectronEnergy_max]))
    wavelengths2_sorted, Ne_final_sorted, pulses1_sorted, Energy_max_sorted = list(itertools.izip(*lists))
    
    if(not EnergyAxis):
        wavelengths2_sorted_show = np.multiply(1e9, wavelengths2_sorted)
    else: #if our wavelength is actually an energy
        wavelengths2_sorted_show = wavelengths2_sorted
    wpi_sorted = np.divide(Ne_final_sorted, np.multiply(1e0, pulses1_sorted)) #NOTE: this is limited to identical pulse durations!
    
    # This loop is run len(Efield2_uniq) times
    #for Efield2_compare in FieldOfInterest: #{{{
    print(Header+"Strange there are 1 FieldOfInterest, but we are getting once more in this loop. Efield2_compare: "+str(Efield2_compare*1E-9))
    #LocalCondition = False
    #if(Efield2_compare==0):
        #LocalCondition = (Efield2_value == Efield2_compare)
    #else:
        #LocalCondition = (abs(Efield2_value - Efield2_compare)/Efield2_value <= FieldTolerance) #Many loops should not induce one new plot each time. 
    #if(LocalCondition): #{{{ #If Efield2_value is close enough to a field we would like to see, then we plot one new line. If we already had the value accepted once, we plot again?? No... So we should construct a set of data by merging all these values. 
    print Header+"Available data for the specific field E2="+str(Efield2_compare*1E-9)+" V/nm."
    tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9 = CountTheNumberOfLaserParameters(PlottingDoublePulses_Efield)
    #### PLOTTING
    
    #if(not EnergyAxis):
    #plot_field, = ax1.plot(wavelengths2_sorted_show, wpi_sorted,'-^', antialiased=antialias, color=Colors[Efield2_compare], linewidth=LineWidthCoeff_abs*LineWidths[Efield2_compare], label=r'$\widetilde{w}_{\mathrm{PI}}$, I='+str(FieldToIntensity(Efield2_value))+r' W/m$^2$') #BACKUP
    #plot_field, = ax1.plot(wavelengths2_sorted_show, wpi_sorted,'-^', antialiased=antialias, label=r'$\widetilde{w}_{\mathrm{PI}}$, I= ? W/m$^2$')
    plot_field, = ax1.plot(wavelengths2_sorted_show, wpi_sorted,'-^', label=r'$\widetilde{w}_{\mathrm{PI}}$, $\lambda_1$='+str(wavelength1*1E9)+' nm, $E_1$='+str(Field1*1E-9)+', $E_2$='+str(Efield2_compare*1E-9)+r' V/nm')
    #linewidth=LineWidthCoeff_abs*LineWidths[Efield2_compare] # BUG 
    plot1 += [plot_field]
    
    #plot_field, = ax1.plot(wavelengths1_Keldysh_Show, np.multiply(
            ##PeakWpi_Gruzdev, np.divide(np.add(n0, np.multiply(-1,PeakWpi_Gruzdev*pulses1_sorted[0])), n0)),
                #'--.'+Colors[Efield2_compare], label=r'$w_{\mathrm{PI}}\frac{n_0-n_{\mathrm{exc}}{n_0}$, I='+str(FieldToIntensity(Efield2_value))+r' W/m$^2$')
    #plot1 += [plot_field]
    
    #}}} #end if
        
    del PlottingDoublePulses_Efield
  #}}} #end for
#}}} #end for
  
  # We separate TDDFT and Keldysh predictions now. 
  #ax2 = pplt.subplot(plotB)
  
  #for Efield2_compare in FieldOfInterest: #{{{
    #print Header+"Building comparison with Keldysh-Gruzdev..."
    #order = 50 
    ##wavelengths1_min = np.min(wavelengths1); wavelengths1_max = np.max(wavelengths1)
    #wavelengths1_min = 200e-9; wavelengths1_max = 2100e-9
    #wavelengths1_Keldysh = np.linspace(wavelengths1_min, wavelengths1_max, 500)
    #PeakWpi_Gruzdev = GenerateKeldyshDatabase(Egap_SI, EffectiveMass[Efield2_compare], wavelengths1_Keldysh, Efield2_compare, order)
    ##PeakWpi_Gulley = GenerateKeldyshGulleyDatabase(Egap_SI, EffectiveMass[Efield2_compare], wavelengths1_Keldysh, Efield2_value, order, 1.)
    ##PeakWpi_Gruzdev_Indirect = GenerateKeldyshDatabase(Egap_SI_indirect, meff, wavelengths1_Keldysh, Efield2_value, order)
    
    #PeakWpi_Tunneling_Direct = np.clip(KeldyshTunnelingLimit(Egap_SI, meff_Tunnel, wavelengths1_Keldysh, Efield2_compare, False), 1E36, 1E46)
    ##PeakWpi_Tunneling_Indirect = np.clip(KeldyshTunnelingLimit(Egap_SI_indirect, meff_Tunnel, wavelengths1_Keldysh, Efield2_value), 1E36, 1E46)
    ##PeakWpi_Gruzdev=
    #PeakWpi_Effective = np.power(np.power(PeakWpi_Gruzdev, -1) + np.power(PeakWpi_Tunneling_Direct,-1),-1)
    ## Dirty: choosing energy or wavelengths. Sorting the results with this choice. 
    #if(EnergyAxis):
        #wavelengths1_Keldysh_Show = np.divide(2.*np.pi*c*hbar/e, wavelengths1_Keldysh) #eV
    #else:
        #wavelengths1_Keldysh_Show = np.multiply(1E9, wavelengths1_Keldysh)
    
    ## resorting the results to plot as function of energy if necessary
    #lists = sorted(itertools.izip(*[wavelengths1_Keldysh_Show, PeakWpi_Gruzdev]))
    #wavelengths1_Keldysh_Show, PeakWpi_Gruzdev_Show = list(itertools.izip(*lists)) 
    
    #plot_field, = ax2.plot(wavelengths1_Keldysh_Show, PeakWpi_Gruzdev_Show, '-', antialiased=antialias, color=Colors[Efield2_compare], linewidth=LineWidthCoeff_rel*LineWidthCoeff_abs*LineWidths[Efield2_compare], label=r'$\widetilde{w}_{\mathrm{PI}}$ [Gruzdev], I='+str(FieldToIntensity(Efield2_value))+r' W/m$^2$')
    
    ##plot_field, = ax1.plot(wavelengths1_Keldysh_Show, PeakWpi_Gruzdev_Show, '-', linewidth=LineWidthCoeff_rel*LineWidthCoeff_abs*LineWidths[Efield2_compare], label=r'$\widetilde{w}_{\mathrm{PI}}$ [Gruzdev], I='+str(FieldToIntensity(Efield2_value))+r' W/m$^2$')
    
    #plot1 += [plot_field]
  ##}}} #end for
    
  #if(EnergyAxis): 
    #ax1.set_xscale('log')
  ax1.set_yscale('log')
  ax1.set_ylabel(r'$\tilde{w}_{PI}$ (m$^{-3}$.s$^{-1}$)')
  #ax2.set_yscale('log')
  
  if(ShowDensities): #{{{
    ax12 = ax1.twinx()
    #pplt.title('Final density as function of laser wavelength')
    print Header+"Total number of datasets: "+str(len(set(Efield2_uniq)))
    
    for Efield2_value in Efield2_uniq: #for each value of Efield1, we make a new search, and plot one new line in the chart
        #PlottingDoublePulses_Efield = FilterDatabase(PlottingDoublePulses, Efield2_value, IndexTable['Efield1']) #select the local Efield
        # TODO: how to merge data with values close to each others? Maybe introduce tolerance of 10% ?
        PlottingDoublePulses_Efield = FilterDatabaseLowerThan(PlottingDoublePulses, Efield2_value*(1.+FieldTolerance), IndexTable['Efield2'])
        PlottingDoublePulses_Efield = FilterDatabaseGreaterThan(PlottingDoublePulses_Efield, Efield2_value*(1.-FieldTolerance), IndexTable['Efield2'])
        
        #print PlottingDoublePulses_Efield
        ## prints the current state of selection and stores the values of sweeped parameters
        runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingDoublePulses_Efield)
        
        if(EnergyAxis): #very dirty!
          wavelengths2=np.divide(2.*np.pi*c*hbar/e, wavelengths2)
        
        lists = sorted(itertools.izip(*[wavelengths2, Ne_final, ElectronEnergy_max])) #NOTE: Ne_max makes more sense than Ne_final
        wavelengths2_sorted, Ne_final_sorted, Energy_max_sorted = list(itertools.izip(*lists))
        if(not EnergyAxis): 
            wavelengths2_sorted = np.multiply(1e9, wavelengths2_sorted)
        
        ax12.set_ylabel(r'Final density $n_{\mathrm{exc}}$ (m$^{-3}$)')
        #print Efield2_value
        # TODO: PUT here a select case which attributes a new color to a new field interval
        # OR USE A DICTIONNARY. 
        
        # SAME CODE, just more clever. 
        for Efield2_compare in FieldOfInterest:
            if(abs(Efield2_value - Efield2_compare)/Efield2_value < FieldTolerance):
                #if(not EnergyAxis):
                plot_field, = ax12.plot(wavelengths2_sorted, Ne_final_sorted, '^', color=Colors[Efield2_compare], label=r'$E_{peak}=$'+str(round(Efield2_compare*1E-9))+' V/nm')

                #plot_field, = ax12.plot(wavelengths2_sorted, Ne_final_sorted, '^', label=r'$E_{peak}=$'+str(round(Efield2_compare*1E-9))+' V/nm')
                plot1 += [plot_field]
                #, I='+str(FieldToIntensity(Efield2_value))+r' W/m$^2$')
                break
        
        #if(1e8*(1.- FieldTolerance) < Efield2_value < 1e8*(1.+ FieldTolerance)):
            ##print Header+"Dataset for the field: "+str(Efield2_value*1e-9)+' V/nm'
            #plot_field, = ax12.plot(wavelengths2_sorted, Ne_final_sorted, '^r', label=r'$E_{peak}=$'+str(round(Efield2_value*1E-9))+' V/nm')
            #plot1 += [plot_field]
            ##, I='+str(FieldToIntensity(Efield2_value))+r' W/m$^2$')
            #break
        #if(1e9*(1.- FieldTolerance) < Efield2_value < 1e9*(1.+ FieldTolerance)):
            ##print Header+"Dataset for the field: "+str(Efield2_value*1e-9)+' V/nm'
            #plot_field, = ax12.plot(wavelengths2_sorted, Ne_final_sorted, '^b', label=r'$E_{peak}=$'+str(round(Efield2_value*1E-9))+' V/nm') 
            #plot1 += [plot_field]
            ##, I='+str(FieldToIntensity(Efield2_value))+r' W/m$^2$')
            #break
        #if(5e9*(1.- FieldTolerance) < Efield2_value < 5e9*(1.+ FieldTolerance)):
            ##print Header+"Dataset for the field: "+str(Efield2_value*1e-9)+' V/nm'
            #plot_field, = ax12.plot(wavelengths2_sorted, Ne_final_sorted, '^g', label=r'$E_{peak}=$'+str(round(Efield2_value*1E-9))+' V/nm')
            #plot1 += [plot_field]
            ##, I='+str(FieldToIntensity(Efield2_value))+r' W/m$^2$')
            #break
        #if(1e10*(1.- FieldTolerance) < Efield2_value < 1e10*(1.+ FieldTolerance)):
            ##print Header+"Dataset for the field: "+str(Efield2_value*1e-9)+' V/nm'
            #plot_field, = ax12.plot(wavelengths2_sorted, Ne_final_sorted, '^k', label=r'$E_{peak}=$'+str(round(Efield2_value*1E-9))+' V/nm')
            #plot1 += [plot_field]
            ##, I='+str(FieldToIntensity(Efield2_value))+r' W/m$^2$')
            #break
        del PlottingDoublePulses_Efield
  #}}}  
  labels1 = [l.get_label() for l in plot1]
  ax1.legend(plot1, labels1, loc='best') #outside right, middle
  #ax1.legend(loc=3, borderaxespad=0., bbox_to_anchor=(1.05, 0.05))
  if(not EnergyAxis): 
    ax1.set_xlabel('Laser wavelength (nm)')
    ax1.set_xlim([200.,2100.])
    #ax2.set_xlabel('Laser wavelength (nm)')
    #ax2.set_xlim([200.,2100.])
  else: 
    ax1.set_xlabel('Photon energy (eV)')
    ax1.set_xlim([0.5,4.0])
    #ax2.set_xlabel('Photon energy (eV)')
    #ax2.set_xlim([0.5,4.0])
  ax1.set_ylim([1E35,1E43])
  #ax2.set_ylim([1E35,1E43])
  #ax1.set_xscale('log')
  #ax1.set_yscale('log')
  #ax1.grid()
  
  ResonantWavelength = h*c/Egap_SI
  print "** Adding vertical bar for the resonant transitions at: "+str(round(1e9*ResonantWavelength))+"nm."
  
  if(not EnergyAxis):
    ax1.axvline(x=1e9*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax1.axvline(x=1e9*2*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax1.axvline(x=1e9*3*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax1.axvline(x=1e9*4*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax1.axvline(x=1e9*wavelength1, color='k', linestyle='--', linewidth=1)
    
    #ax2.axvline(x=1e9*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    #ax2.axvline(x=1e9*2*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    #ax2.axvline(x=1e9*3*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    #ax2.axvline(x=1e9*4*ResonantWavelength, color='k', linestyle='--', linewidth=1)
  else:
    ax1.axvline(x=Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax1.axvline(x=1./2.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax1.axvline(x=1./3.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax1.axvline(x=1./4.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax1.axvline(x=2.*np.pi*c/wavelength1/e, color='r', linestyle='-', linewidth=2)
    
    #ax2.axvline(x=Egap_SI/e, color='k', linestyle='--', linewidth=1)
    #ax2.axvline(x=1./2.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    #ax2.axvline(x=1./3.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    #ax2.axvline(x=1./4.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
  
  ## Saving this first double figure
  #pplt.tight_layout()
  filename = 'Metaplot-DoublePulse-BiColor-EffectOfWavelength'+str(round(wavelength1*1E9,0))+'nm-Efield1_Tol-'+str(round(Field1*1E-9,1))+'Vnm-Keldysh-'+str(int(pulseDuration*1E15))+"-CEP"+str(round(CEP2_request*180./np.pi,0))
  pplt.savefig(filename+'.eps')
  pplt.savefig(filename+'.png')
  
  ## Imports the optical data from Palik database
  # Returns the spectrum (in Re(eps), Im(eps)) of a material on the given grid of wavelengths
  def ImportPalikDatabase_epsilon_fromNK(wavelengths, DataFile="Si-Palik"): #{{{
      DataArray           = np.loadtxt(DataFile, delimiter="\t", skiprows=4)
      wavelengths_source  = DataArray[:,0]
      if(np.min(np.diff(wavelengths_source))==0e0): 
          print Header+"Error in wavelengths_source. Careful! Palik book combines different optical data, hence there can be overlap, sometimes! It must be removed to be imported properly here. " 
          #print np.diff(wavelengths_source)
          exit()
      OpticalN            = DataArray[:,1]
      OpticalK            = DataArray[:,2]
      order = 1
      fni = InterpolatedUnivariateSpline(wavelengths_source*1E-10, OpticalN, k=order)
      fki = InterpolatedUnivariateSpline(wavelengths_source*1E-10, OpticalK, k=order)
      OpticalN_new = fni(wavelengths)
      OpticalK_new = fki(wavelengths)
      epsilon = np.power(OpticalN_new+np.multiply(1j, OpticalK_new), 2)
      #epsR = np.real(epsilon); epsC = np.imag(epsilon)
      return epsilon
  #}}}
  
  ## Computes the free carrier energy contained in the volume <SampleVolume>.  
  def EnergyFromFreeCarrierAbsorption(pulseDuration, wavelengths1_sorted, Efield, Ne_final_sorted, w_PI_sorted, Egap_SI = 2.56*e, meff = 0.2226, SampleVolume=2 * 5E28**-1): #{{{
    # dEnergy/dt = alpha * intensity + ( n har omega - Egap ) x sigma_n I^n
    # where we approach sigma_n I^n by pulse-averaged w_PI 
    # Energy ~ alpha * intensity + ( n har omega - Egap ) x w_PI * tau
    omega   = np.divide(2. * np.pi * c, wavelengths1_sorted)
    #print Header+"omega: "+str(omega) #Valid
    nu = 6e-15**-1 #1.1 fs: Solowski-Tinten adjustment. Can also use the one from Shunsuke Sato, PRB 2015 (Numerical pump probe)
    print "nu: "+str(nu) #Valid
    #epsilon_Si0 = 13.64 + 0.048j #TODO: get the optical data in extended spp theory for Si
    epsilon_Si0 = ImportPalikDatabase_epsilon_fromNK(wavelengths1_sorted) ##TODO: replace by LDA !
    print "Dielectric permittivity of Si at 800 nm : "+str(ImportPalikDatabase_epsilon_fromNK(800e-9))
    omega_p2    = np.multiply(Ne_final_sorted, e**2 / m_e / meff / epsilon_0)
    #print "omega_p: "+str(np.sqrt(omega_p2)) # Valid
    # NOTE: THIS IS A CHOICE. WE JUST HAVE TO BE CONSISTENT ON WHICH PROCESS WE WISH FOR. 
    Epsilon = 1E0 - np.divide(omega_p2, np.power(omega, 2)) * np.divide(1.,1. + np.divide(1.j * nu, omega))
    #Epsilon = epsilon_Si0 - np.divide(omega_p2, np.power(omega, 2)) * np.divide(1.,1. + np.divide(1.j * nu, omega))
    #NumberOfPhotons = Egap_SI / (hbar * omega)
    NumberOfPhotons = np.divide(Egap_SI/hbar, omega)
    print NumberOfPhotons
    #NumberOfPhotons_round = np.round(NumberOfPhotons,0)
    NumberOfPhotons_round = np.ceil(NumberOfPhotons)
    print Header+"sqrt(Epsilon): "+str(np.min(np.sqrt(Epsilon)))+", "+str(np.max(np.sqrt(Epsilon))) #Valid
    alpha                = np.multiply(2.* omega/c, np.imag(np.sqrt(Epsilon))) #valid
    EnergyExcess         = np.multiply(NumberOfPhotons_round * hbar * omega, Ne_final_sorted) * SampleVolume / e  #result in eV
    EnergyBremsstrahlung = np.multiply( omega * epsilon_0 * (Epsilon.imag), np.multiply(Efield, Efield)) * SampleVolume / e  * pulseDuration
    EnergyBremsstrahlung_approx =  EnergyBremsstrahlung #NOTE: not used for now. Another collision freq? 
    #EnergyBremsstrahlung_approx = np.multiply( 4.*np.pi / c * np.divide(Epsilon.imag, (np.sqrt(Epsilon)).real), )) ) * SampleVolume / e
    print Header+"FcrEnergy (eV): "+str(np.min(EnergyBremsstrahlung))+", "+str(np.max(EnergyBremsstrahlung))
    print Header+"Excess energy (n hbar omega) (eV): "+str(np.min(EnergyExcess))+", "+str(np.max(EnergyExcess))
    EnergyFreeCarrier    = EnergyExcess + EnergyBremsstrahlung
    EnergyFreeCarrier_approx = EnergyExcess + EnergyBremsstrahlung_approx
    return EnergyFreeCarrier, EnergyFreeCarrier_approx
  #}}}
  
  ## TODO: TDDFT vs Drude should become another figure. We don't need to have them side by side. 
  del labels1; del plot1
  plt.figure(figsize=(12,6))
  
  plot1=[]
  
  ax3 = plt.subplot(111)
  #plt.title('Irreversibility factor')
  #pplt.plot(FieldToIntensity(Efield1), irreversibility, 'o', label='Irreversibility')
  for Efield2_compare in FieldOfInterest: #for each value of Efield1, we make a new search, and plot one new line in the chart
    #PlottingDoublePulses_Efield = FilterDatabase(PlottingDoublePulses, Efield2_value, IndexTable['Efield1']) #select the local Efield
    PlottingDoublePulses_Efield = FilterDatabaseLowerThan(PlottingDoublePulses, Efield2_compare*(1.+FieldTolerance), IndexTable['Efield2'])
    PlottingDoublePulses_Efield = FilterDatabaseGreaterThan(PlottingDoublePulses_Efield, Efield2_compare*(1.-FieldTolerance), IndexTable['Efield2'])
    #print Header+"Dataset for the field: "+str(Efield2_value*1e-9)+' V/nm'
    #print PlottingDoublePulses_Efield
    # TODO: as fields are all slighly different, we should project on the neareast neighboring value, up to 1E8 V/m. 
    
    ## prints the current state of selection and stores the values of sweeped parameters
    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingDoublePulses_Efield)
    
    # 
    if(EnergyAxis):
        wavelengths2_Show = np.divide(2.*np.pi*c*hbar/e, wavelengths2) #eV
    else: 
        wavelengths2_Show = np.multiply(1e9, wavelengths2)
    
    lists = sorted(itertools.izip(*[wavelengths2_Show, wavelengths2, Ne_final, W_pi, ElectronEnergy_max])) #NOTE: Ne_max makes more sense than Ne_final
    wavelengths2_sorted_show, wavelengths2_sorted, Ne_final_sorted, W_pi_sorted, Energy_max_sorted = list(itertools.izip(*lists))
    #print Efield2_value
    
    #for Efield2_compare in FieldOfInterest: 
    #LocalCondition = False
    #if(Efield2_compare==0):
        #LocalCondition = (Efield2_value == Efield2_compare)
    #else:
        #LocalCondition = (abs(Efield2_value - Efield2_compare)/Efield2_compare <= FieldTolerance)
    #if(LocalCondition):  #if we have captured the fields we are interested in, we plot  
    #ax3.plot(wavelengths2_sorted_show, Energy_max_sorted, '-', linewidth=LineWidthCoeff_abs*LineWidths[Efield2_compare], color=Colors[Efield2_compare])
    plot_field, = ax3.plot(wavelengths2_sorted_show, Energy_max_sorted, '-^', label=r'$\lambda_1=$'+str(wavelength1*1E9)+' nm, $E_1=$'+str(Field1*1E-9)+r' V/nm, $E_2=$'+str(Efield2_compare*1E-9)+" V/nm.")
    print(Header+"Plotting energy: Adding one plot...")
    print(wavelengths2_sorted_show)
    print(Energy_max_sorted)
    plot1 += [plot_field]
    #ax3.plot(wavelengths2_sorted_show, Energy_max_sorted, '-')
            
    del PlottingDoublePulses_Efield
    
  #pplt.xlabel(r'Intensity (W/m$^2$)')
  if(not EnergyAxis):
    ax3.set_xlabel('Laser wavelength (nm)')
  else: 
    ax3.set_xlabel('Photon energy (eV)')
  #pplt.ylabel('Irreversibility factor')
  ax3.set_ylabel('Electron excess energy (eV)')
  #pplt.ylim((0,1))
  #ax2.set_xscale('log')
  
  ## TODO: this can be removed, since irreversibility is gauge-dependent. 
  #if(PlotIrreversibility): #{{{
    #ax22 = ax2.twinx()
    #for Efield2_value in Efield2_uniq: #for each value of Efield1, we make a new search, and plot one new line in the chart
        ##PlottingDoublePulses_Efield = FilterDatabase(PlottingDoublePulses, Efield2_value, IndexTable['Efield1']) #select the local Efield
        #PlottingDoublePulses_Efield = FilterDatabaseLowerThan(PlottingDoublePulses, Efield2_value*(1.+FieldTolerance), IndexTable['Efield1'])
        #PlottingDoublePulses_Efield = FilterDatabaseGreaterThan(PlottingDoublePulses_Efield, Efield2_value*(1.-FieldTolerance), IndexTable['Efield1'])
        ##print Header+"Dataset for the field: "+str(Efield2_value*1e-9)+' V/nm'
        ##print PlottingDoublePulses_Efield
        ## TODO: as fields are all slighly different, we should project on the neareast neighboring value, up to 1E8 V/m. 
        
        ### prints the current state of selection and stores the values of sweeped parameters
        #runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingDoublePulses_Efield)
        
        #lists = sorted(itertools.izip(*[wavelengths1, Ne_final, Energy_max, irreversibility])) #NOTE: Ne_max makes more sense than Ne_final
        #wavelengths1_sorted, Ne_final_sorted, Energy_max_sorted, irreversibility_sorted = list(itertools.izip(*lists))
        #wavelengths1_sorted_show = np.multiply(1e9, wavelengths1_sorted)
        
        #for Efield_compare in FieldOfInterest:
            #if(abs(Efield2_value - Efield2_compare)/Efield2_compare < FieldTolerance):
                #ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+'+Colors[Efield2_compare], label='Irreversibility')
        
        ###print Efield2_value
        ### TODO: PUT here a select case which attributes a new color to a new field interval
        ###if(1e8*(1.-FieldTolerance) < Efield2_value < 1e8*(1.+FieldTolerance)):
            ###ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+r', label='Irreversibility')
            ###ax22.plot(wavelengths1_sorted_show, gammaKeldysh(Egap_SI, 0.228, Efield2_value, wavelengths1_sorted), '--r', label='Adiabaticity')
        ##if(1e9*(1.-FieldTolerance) < Efield2_value < 1e9*(1.+FieldTolerance)):
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+b', label='Irreversibility')
            ###ax22.plot(wavelengths1_sorted_show, gammaKeldysh(Egap_SI, 0.228, Efield2_value, wavelengths1_sorted), '--b', label='Adiabaticity')
        ##if(5e9*(1.-FieldTolerance) < Efield2_value < 5e9*(1.+FieldTolerance)):
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+g', label='Irreversibility')
            ###ax22.plot(wavelengths1_sorted_show, gammaKeldysh(Egap_SI, 0.228, Efield2_value, wavelengths1_sorted), '--g', label='Adiabaticity')
        ##if(1e10*(1.-FieldTolerance) < Efield2_value < 1e10*(1.+FieldTolerance)):
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+k', label='Irreversibility')
            ###ax22.plot(wavelengths1_sorted_show, gammaKeldysh(Egap_SI, 0.228, Efield2_value, wavelengths1_sorted), '--k', label='Adiabaticity')
        #del PlottingDoublePulses_Efield
    
    
    #pplt.ylabel(r'Irreversibility $\zeta$')
    #ax22.set_yscale('log')
  ##}}}
  
  # This part adds an analytical estimation based on Drude(n_exc_TDDFT)
  if(PlotFreeCarrierAbsorption):
      print Header+"Now plotting the free carrier energy absorption..."
      for Efield2_value in Efield2_uniq: #for each value of Efield1, we make a new search, and plot one new line in the chart
        #PlottingDoublePulses_Efield = FilterDatabase(PlottingDoublePulses, Efield2_value, IndexTable['Efield1']) #select the local Efield
        PlottingDoublePulses_Efield = FilterDatabaseLowerThan(PlottingDoublePulses, Efield2_value*(1.+FieldTolerance), IndexTable['Efield2'])
        PlottingDoublePulses_Efield = FilterDatabaseGreaterThan(PlottingDoublePulses_Efield, Efield2_value*(1.-FieldTolerance), IndexTable['Efield2'])
        #print Header+"Dataset for the field: "+str(Efield2_value*1e-9)+' V/nm'
        #print PlottingDoublePulses_Efield
        # TODO: as fields are all slighly different, we should project on the neareast neighboring value, up to 1E8 V/m. 
        
        ## prints the current state of selection and stores the values of sweeped parameters
        runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingDoublePulses_Efield)
        
        if(EnergyAxis):
            wavelengths2_Show = np.divide(2.*np.pi*c*hbar/e, wavelengths2) #eV
        else: 
            wavelengths2_Show = np.multiply(1e9, wavelengths2)
    
        lists = sorted(itertools.izip(*[wavelengths2_Show, wavelengths2, Ne_final, W_pi, ElectronEnergy_max])) #NOTE: Ne_max makes less sense than Ne_final
        wavelengths2_sorted_show, wavelengths2_sorted, Ne_final_sorted, W_pi_sorted, Energy_max_sorted = list(itertools.izip(*lists))
    #print Efield2_value
    
        
        #lists = sorted(itertools.izip(*[wavelengths2, Ne_final, W_pi, Energy_max, irreversibility])) #NOTE: Ne_max makes less sense than Ne_final
        #wavelengths2_sorted, Ne_final_sorted, W_pi_sorted, Energy_max_sorted, irreversibility_sorted = list(itertools.izip(*lists))
        #wavelengths2_sorted_show = np.multiply(1e9, wavelengths2_sorted)
        
        for Efield2_compare in FieldOfInterest:
            LocalCondition = False
            if(Efield2_compare==0):
                LocalCondition = (Efield2_value == Efield2_compare)
            else:
                LocalCondition = (abs(Efield2_value - Efield2_compare)/Efield2_compare <= FieldTolerance)
            if(LocalCondition):
                EnergyAnalytical, EnergyAnalytical_approx = EnergyFromFreeCarrierAbsorption(pulses1, wavelengths2_sorted, Efield2_value, Ne_final_sorted, W_pi_sorted, Egap_SI, EffectiveMass[Efield2_compare])
                                
                plot_field, = ax3.plot(wavelengths2_sorted_show, EnergyAnalytical, '--', color=Colors[Efield2_compare], label=r'$E_1=$'+str(Field1*1E-9)+r' V/nm, $E_2=$'+str(Efield2_compare*1E-9)+" V/nm.")
                #ax2.plot(wavelengths2_sorted_show, EnergyAnalytical, '--', label='Free carrier absorption')
                plot1 += [plot_field]
            
        ##if(1e8*(1.-FieldTolerance) < Efield2_value < 1e8*(1.+FieldTolerance)):
            ##print Header+"Found field 0.1 V/nm"
            ##EnergyAnalytical, EnergyAnalytical_approx = EnergyFromFreeCarrierAbsorption(pulses1, wavelengths1_sorted, Efield2_value, Ne_final_sorted, W_pi_sorted, Egap_SI)
            ##ax2.plot(wavelengths1_sorted_show, EnergyAnalytical, '--+r', label='Free carrier absorption')
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+r', label='Irreversibility')
            ##ax22.plot(wavelengths1_sorted_show, gammaKeldysh(Egap_SI, 0.228, Efield2_value, wavelengths1_sorted), '--r', label='Adiabaticity')
        #if(1e9*(1.-FieldTolerance) < Efield2_value < 1e9*(1.+FieldTolerance)):
            #print Header+"Found field 1.0 V/nm"
            #EnergyAnalytical, EnergyAnalytical_approx = EnergyFromFreeCarrierAbsorption(pulses1, wavelengths1_sorted, Efield2_value, Ne_final_sorted, W_pi_sorted, Egap_SI)
            #ax2.plot(wavelengths1_sorted_show, EnergyAnalytical, '--+b', label='Free carrier absorption')
            ##ax2.plot(wavelengths1_sorted_show, EnergyAnalytical_approx, '--b', label='Free carrier absorption')
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+b', label='Irreversibility')
            ##ax22.plot(wavelengths1_sorted_show, gammaKeldysh(Egap_SI, 0.228, Efield2_value, wavelengths1_sorted), '--b', label='Adiabaticity')
        #if(5e9*(1.-FieldTolerance) < Efield2_value < 5e9*(1.+FieldTolerance)):
            #print Header+"Found field 5 V/nm"
            #EnergyAnalytical, EnergyAnalytical_approx = EnergyFromFreeCarrierAbsorption(pulses1, wavelengths1_sorted, Efield2_value, Ne_final_sorted, W_pi_sorted, Egap_SI)
            #ax2.plot(wavelengths1_sorted_show, EnergyAnalytical, '--+g', label='Free carrier absorption')
            ##ax2.plot(wavelengths1_sorted_show, EnergyAnalytical_approx, '--g', label='Free carrier absorption')
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+g', label='Irreversibility')
            ##ax22.plot(wavelengths1_sorted_show, gammaKeldysh(Egap_SI, 0.228, Efield2_value, wavelengths1_sorted), '--g', label='Adiabaticity')
        #if(1e10*(1.-FieldTolerance) < Efield2_value < 1e10*(1.+FieldTolerance)):
            #print Header+"Found field 10 V/nm"
            #EnergyAnalytical, EnergyAnalytical_approx = EnergyFromFreeCarrierAbsorption(pulses1, wavelengths1_sorted, Efield2_value, Ne_final_sorted, W_pi_sorted, Egap_SI)
            #ax2.plot(wavelengths1_sorted_show, EnergyAnalytical, '--+k', label='Free carrier absorption')
            ##ax2.plot(wavelengths1_sorted_show, EnergyAnalytical_approx, '--k', label='Free carrier absorption')
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+k', label='Irreversibility')
            ##ax22.plot(wavelengths1_sorted_show, gammaKeldysh(Egap_SI, 0.228, Efield2_value, wavelengths1_sorted), '--k', label='Adiabaticity')
        del PlottingDoublePulses_Efield
    
    #pplt.ylabel(r'Irreversibility $\zeta$')
    #ax22.set_yscale('log')
      #ax22 = ax2.twinx()
  #pplt.ylabel(r'Adiabaticity $\gamma$, Irreversibility $\zeta$')
  #ax22.set_ylim((0,1))
  labels1 = [l.get_label() for l in plot1]
  ax3.legend(plot1, labels1, loc='best') #outside right, middle
  #ax3.set_xscale("log")
  ax3.set_yscale("symlog", linthreshy=1e-9)
  #ax2.set_ylim((1E-9,5E2))
  #ax22.grid()
  if(not EnergyAxis): 
    pplt.xlim((200.,2100.))
  else: 
    pplt.xlim((0.5,7.))
  
  
  #ResonantWavelength = h*c/Egap_SI
  print "** Adding vertical bar for the resonant transitions at: "+str(round(1e9*ResonantWavelength))+"nm."
  
  if(not EnergyAxis):
    ax3.axvline(x=1e9*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax3.axvline(x=1e9*2*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax3.axvline(x=1e9*3*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax3.axvline(x=1e9*4*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax3.axvline(x=1E9*wavelength1, color='r', linestyle='-', linewidth=2)
  else:
    ax3.axvline(x=Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax3.axvline(x=1./2.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax3.axvline(x=1./3.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax3.axvline(x=1./4.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax3.axvline(x=2.*np.pi*c/wavelength1/e, color='r', linestyle='-', linewidth=2)
  #ax2.legend(loc=2, borderaxespad=0.)
  #ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
  #ax22.legend(loc=1, borderaxespad=0.)

  # Saving the plot
  #pplt.tight_layout()
  #filename = 'Metaplot-DoublePulse-BiColor-EffectOfWavelength'+str(round(wavelength1*1E9,0))+'nm-Efield1_Tol-'+str(round(Efield1*1E-9,0))+'Vnm-Keldysh-'+str(int(pulseDuration*1E15))
  filename='Metaplot-DoublePulse-BiColor-EffectOfWavelength'+str(round(wavelength1*1E9,0))+'nm-Efield1_Tol-'+str(round(Field1*1E-9,1))+'Vnm-Drude-'+str(int(pulseDuration*1E15))
  pplt.savefig(filename+'.eps')
  pplt.savefig(filename+'.png')
  print(Header+"Bicolor_plotEffectOfWavelength_Tolerance: successfully wrote the bicolor Stark analysis.")
  if(ShowFinalPlot==1):
    pplt.show()

  del ax1, PlottingDoublePulses
  #}}}

## Plots Nexc | Energy (CEP) for bicolor simulations. Selection is made with function parameters. A 1 percent tolerance is used on the <Efield1> and <Efield2> amplitudes. The rest is exactly tested. 
def Bicolor_plotEffectOfCEP(PlottingTable, Color1=800e-9, duration1=20e-15, FieldAmplitude1=1E9, Color2=1600e-9, duration2=20e-15, FieldAmplitude2=1E9, PulseShape=3, pulseRiseTime=0, ShowFinalPlot=0, Egap_SI=2.56*e):
  print "#====== Plot: DOUBLE COLOR PULSES, EFFECT OF THE CEP, SAME POLARIZATIONS, SAME AMPLITUDE: final density as function of field ======"
  #Filter: double pulses, same wavelengths, same polarization
  
  # Simple pulse, single color
  PlottingSinglePulses = PlottingTable #get the full database
  PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
  PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select double pulses.
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, Color1, IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts. #TODO: this field should be making a loop in available wavelengths given in wavelengths1_uniq.
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, Color2, IndexTable['wavelengths2']) #keeps only 800 nm wavelenghts. #TODO: this field should be making a loop in available wavelengths given in wavelengths1_uniq.
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0, IndexTable['FreezeHXC']) #keeps only non-frozen XC. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, PulseShape, IndexTable['PulseType']) #keeps only non-frozen XC. 
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, duration1, IndexTable['pulses1']) #removes fields where tau != pulseDuration
  PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, duration2, IndexTable['pulses2']) #removes fields where tau != pulseDuration
  PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, FieldAmplitude1*1.01, IndexTable['Efield1'])
  PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, FieldAmplitude1*0.99, IndexTable['Efield1'])
  PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, FieldAmplitude2*1.01, IndexTable['Efield2'])
  PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, FieldAmplitude2*0.99, IndexTable['Efield2'])
  
  if(pulseRiseTime != 0): 
    PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, pulseRiseTime, IndexTable['RiseTime']) #removes fields where tau != pulseDuration

  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = CountTheNumberOfLaserParameters(PlottingSinglePulses)
  
  pplt.figure(figsize=(8,8./GoldNumber))
  ax1 = pplt.subplot(111) #111 or plotA if using plotB
  pplt.title(r'$n_{\mathrm{exc}}(\phi_{\mathrm{CEP}})$')
  # Number of sweeped parameters: 3 (Efield, Delay, CrossPolarization)
  #pplt.scatter(FieldToIntensity(Efield1), Ne_final, c=Delay, s=(30.+CrossPolarization*50.), label=r'$N_{exc}$[Field]')
  plot1 = [] #np.empty((0,0))
  #for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
  #PlottingSinglePulses_CEP = FilterDatabase(PlottingSinglePulses, CEP1_value, IndexTable['CEP1']) #select the local CEP
  
  ## prints the current state of selection and stores the values of sweeped parameters
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses)
  
  #plot_field, = ax1.plot(FieldToIntensity(Efield1), Ne_final, '-', label=r'$N_{exc}^{final}$, $\phi_{CE}$='+str(round(CEP1_value*180./pi)))
  #lists = sorted(itertools.izip(*[FieldToIntensity(Efield1), Ne_final]))
  lists = sorted(itertools.izip(*[CEP2, Ne_final, W_pi, ElectronEnergy_max, Efield2]))
  CEP2_sorted, Ne_final_sorted, W_pi_sorted, Energy_max_sorted, Efield2_sorted = list(itertools.izip(*lists))
  plot_field, = ax1.plot(np.multiply(CEP2_sorted,180./pi), Ne_final_sorted, 'r-+', label=r'$n_{\mathrm{exc}}^{\mathrm{final}}$, $E^{\mathrm{peak}}_{1,2}=$'+str(round(FieldAmplitude1*1E-9))+' V/nm.')
  print plot_field
  print Header+"** Debug info:"+str(np.shape(plot_field))
  plot1 += [plot_field] # does cumulate each CEP now. Type: Line2D, not np.array. 
    
  #ax1.legend(loc=4, borderaxespad=0.) #bbox_to_anchor=(1.05, 1), 
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
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses)
  
  #plot_field, = ax1.plot(FieldToIntensity(Efield1), Ne_final, '-', label=r'$N_{exc}^{final}$, $\phi_{CE}$='+str(round(CEP1_value*180./pi)))
  #lists = sorted(itertools.izip(*[FieldToIntensity(Efield1), Ne_final]))
  lists = sorted(itertools.izip(*[CEP2, Ne_final, W_pi, ElectronEnergy_max, pulses2, Efield2]))
  CEP2_sorted, Ne_final_sorted, W_pi_sorted, Energy_max_sorted, pulses1_sorted, Efield2_sorted  = list(itertools.izip(*lists))
  #plt.plot(np.multiply(CEP2_sorted,180./pi), np.divide(Ne_final_sorted, pulses1_sorted), '-') #, label=r'$N_{exc}/\tau$') #, I='+str(FieldToIntensity(Efield1_value))+r' W/m$^2$')
  plot_field, = ax12.plot(np.multiply(CEP2_sorted,180./pi), Energy_max_sorted, 'b-', label=r'$\xi_{\mathrm{el}}^{\mathrm{TDDFT}}$ , $E^{\mathrm{peak}}_{1,2}=$'+str(round(FieldAmplitude2*1E-9,0))+' V/nm')
  plot1 += [plot_field]
  
  #ax12.set_yscale('log')
  ax12.set_ylabel(r'$\xi_{\mathrm{el}}^{\mathrm{TDDFT}}$ (eV)')
  #ax.ylabel('Virtual exc. factor')
  
  labels1 = [l.get_label() for l in plot1]
  ax1.legend(plot1, labels1, loc='best', fancybox=True, framealpha=0)
  
  
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
  ExportName = "EffectOfCEP-E1E2-"+str(round(FieldAmplitude1*1E-9,0))+'-tau1tau2-'+str(round(duration1*1E15,0))
  plt.savefig('Metaplot-Bicolor-'+ExportName+'.eps')
  plt.savefig('Metaplot-Bicolor-'+ExportName+'.png')
  if(ShowFinalPlot==1):
    pplt.show()

  del ax1, PlottingSinglePulses, wavelength_uniq
#}}}

## Plots the N_exc(wavelength, Efield) (2D plot)
# Not usable, unless you have computed a continuum of wavelengths AND intensities (crazy amount of computer time - 1M core-hour per figure)
# VALID, but of course interference of 2 waves is oscillating very much. 
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
      runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingSinglePulses_Efield_wavelength)
      
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


# Plots (Ne_final | Wpi | Energy | E**2) (Efield1, Efield2) for pulses of fixed wavelengths, same pulse duration, same CEP. 
# THIS FIGURE IS THE ONE THAT CAN BE COMPARED WITH VLADIMIR ZHUKOV. 
# From this function, we perform a  wish to GENERATE a 1D section in this figure to analyze the role of Stark effect in bicolor absorption.
def plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1=800e-9, wavelength2=1600e-9, PulseDuration=20e-15, PulseType=3, RiseCycles=2, WhatToPlot=0, Egap_SI=2.56*e, RemoveTunneling=True): #{{{
  LogWpi=False                                                                                                                                                                                                                
  LogPlotXY=False
  gammaLimit = 1E-10
  meff = 0.2226 #effective mass of electron [Waroquiers, D. et al, Physical Review B, 2013, 87, 075121]
  CEP2_request=-pi/2.
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
  
  if(CEP2_request >= 0):
      PlottingDoublePulses = FilterDatabaseLowerThan(PlottingDoublePulses,   (CEP2_request)*(1.1), IndexTable['CEP2'])
      PlottingDoublePulses = FilterDatabaseGreaterThan(PlottingDoublePulses, (CEP2_request)*(0.9), IndexTable['CEP2']) #keeps only CEP = 0. 
  else:
      PlottingDoublePulses = FilterDatabaseGreaterThan(PlottingDoublePulses, (CEP2_request)*(1.1), IndexTable['CEP2'])
      PlottingDoublePulses = FilterDatabaseLowerThan(PlottingDoublePulses,   (CEP2_request)*(0.9), IndexTable['CEP2']) #keeps only CEP = 0. 
  # TODO: we shall clean CEP2 as well, as it must influence the results. 
  
  #PlottingDoublePulses = FilterDatabase(PlottingDoublePulses, CEP2, IndexTable['CEP2']) # Must have same CEP
  
  if(RiseCycles != 0): #0 selects all the possible types
    PlottingDoublePulses = FilterDatabase(PlottingDoublePulses, RiseCycles, IndexTable['RiseTime']) # Selecting Rise Time
  #### 2D PLOTTING Density(Efield1, Efield2) for one given set of laser parameters. For each combination, we produce one plot like this. 
  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq, wavelength2_uniq, CEP2_uniq, Efield2_uniq, pulses2_uniq= CountTheNumberOfLaserParameters(PlottingDoublePulses)
  print Header+"Total number of datasets: "+str(len(set(Efield1_uniq))*len(set(Efield2_uniq)))
  
  TotalField2D      = np.empty((0,3)) #preparing the 2D plot
  DensityMax2D      = np.empty((0,3)) #preparing the 2D plot
  DensityMax2D_grad = np.empty((0,3)) #preparing the 2D plot
  CurrentMax2D      = np.empty((0,3)) #preparing the 2D plot
  EnergyMax2D       = np.empty((0,3)) #preparing the 2D plot
  WpiMax2D          = np.empty((0,3)) #preparing the 2D plot
  
  for Efield1_i in Efield1_uniq: #{{{ #selecting amplitude1
    PlottingDoublePulses_SelectedField1 = FilterDatabase(PlottingDoublePulses, Efield1_i, IndexTable['Efield1'])
    for Efield2_i in Efield2_uniq: #{{{ #selecting amplitude2
      # Selecting Field2
      PlottingDoublePulses_SelectedField2 = FilterDatabase(PlottingDoublePulses_SelectedField1, Efield2_i, IndexTable['Efield2'])
      
      # Constructing the density table
      runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingDoublePulses_SelectedField2)
      print Header+"Field1: "+str(Efield1)+" m, Field2: "+str(Efield2)+" m."
      
      if(len(Ne_final) >= 1): 
        Field1_gamma = gammaKeldysh(Egap_SI, meff, Efield1[0], wavelengths1[0])
        Field2_gamma = gammaKeldysh(Egap_SI, meff, Efield2[0], wavelengths2[0])
      
        # Tests if regime is tunneling
        if(RemoveTunneling): 
            Field1_Multiphoton = (Field1_gamma > gammaLimit)
            Field2_Multiphoton = (Field2_gamma > gammaLimit)
        else: #we dont check
            Field1_Multiphoton = True
            Field2_Multiphoton = True
      
        if(Field1_Multiphoton and Field2_Multiphoton): #We filter single-pulses in tunneling regime
            ToBeAddedField  =np.array([Efield1[0], Efield2[0], np.power(Efield1[0]+Efield2[0],2)]) #builds the total peak squared field
            ToBeAdded       =np.array([Efield1[0], Efield2[0], Ne_final[0]])                       #BUG: removes some results, not clean enough!
            ToBeAddedGrad   =np.array([.5*c*epsilon_0*np.power(Efield1[0],2), .5*c*epsilon_0*np.power(Efield2[0],2), Ne_final[0]]) #meant to find d(Nexc)/dI_1,2
            ToBeAddedCurrent=np.array([Efield1[0], Efield2[0], TotalCurrent_final[0]]) #BUG: removes some results, not clean enough!
            ToBeAddedEnergy =np.array([Efield1[0], Efield2[0], ElectronEnergy_max[0]]) #BUG: removes some results, not clean enough!
            if(LogWpi):
                ToBeAddedWpi    =np.array([Efield1[0], Efield2[0], np.log10(W_pi[0])]) #BUG: removes some results, not clean enough!
            else: 
                ToBeAddedWpi    =np.array([Efield1[0], Efield2[0], W_pi[0]]) #BUG: removes some results, not clean enough!
            #print ToBeAdded
            DensityMax2D_grad=np.vstack((DensityMax2D_grad, ToBeAddedGrad)) 
            TotalField2D     =np.vstack((TotalField2D     , ToBeAddedField))
            DensityMax2D     =np.vstack((DensityMax2D     , ToBeAdded))
            CurrentMax2D     =np.vstack((CurrentMax2D     , ToBeAddedCurrent))
            EnergyMax2D      =np.vstack((EnergyMax2D      , ToBeAddedEnergy))
            WpiMax2D         =np.vstack((WpiMax2D         , ToBeAddedWpi))
            print "Cleaning..."
            del ToBeAdded, ToBeAddedCurrent, ToBeAddedEnergy, ToBeAddedWpi, ToBeAddedField
    #}}}
  #}}}
  PulseTypeName = int(np.asfarray(PulseType)[0])
  # Stack is built. Export, show and plot on a heat map. 
  print "** Info: exporting the sub-plot2D: Ne_final(Field1, Field2)..."
  print "** Info: useful to visualize available data with GNUplot. "
  filename_field   = "Bicolor-"+str(PulseTypeName)+"-Color1-"+str(round(wavelength1*1E9))+"-Color2-"+str(round(wavelength2*1E9))+"-"+str(int(PulseDuration*1e15))+"fs-EPeakWithE1E2"
  filename         = "Bicolor-"+str(PulseTypeName)+"-Color1-"+str(round(wavelength1*1E9))+"-Color2-"+str(round(wavelength2*1E9))+"-"+str(int(PulseDuration*1e15))+"fs-NeFinalWithE1E2"
  filename_gradNexc= "Bicolor-"+str(PulseTypeName)+"-Color1-"+str(round(wavelength1*1E9))+"-Color2-"+str(round(wavelength2*1E9))+"-"+str(int(PulseDuration*1e15))+"fs-NexcWithI1I2"
  filename_current = "Bicolor-"+str(PulseTypeName)+"-Color1-"+str(round(wavelength1*1E9))+"-Color2-"+str(round(wavelength2*1E9))+"-"+str(int(PulseDuration*1e15))+"fs-CurrentFinalWithE1E2"
  filename_energy  = "Bicolor-"+str(PulseTypeName)+"-Color1-"+str(round(wavelength1*1E9))+"-Color2-"+str(round(wavelength2*1E9))+"-"+str(int(PulseDuration*1e15))+"fs-EnergyFinalWithE1E2"
  filename_wpi     = "Bicolor-"+str(PulseTypeName)+"-Color1-"+str(round(wavelength1*1E9))+"-Color2-"+str(round(wavelength2*1E9))+"-"+str(int(PulseDuration*1e15))+"fs-WpiWithE1E2"
  # Sorting the data before storing and plotting
  TotalField2Dsorted = TotalField2D[np.lexsort((TotalField2D[:,2], TotalField2D[:,1], TotalField2D[:,0]))]
  DensityMax2Dsorted = DensityMax2D[np.lexsort((DensityMax2D[:,2], DensityMax2D[:,1], DensityMax2D[:,0]))]
  DensityMax2DgradSorted = DensityMax2D_grad[np.lexsort((DensityMax2D_grad[:,2], DensityMax2D_grad[:,1], DensityMax2D_grad[:,0]))]
  CurrentMax2Dsorted = CurrentMax2D[np.lexsort((CurrentMax2D[:,2], CurrentMax2D[:,1], CurrentMax2D[:,0]))]
  EnergyMax2Dsorted  = EnergyMax2D[np.lexsort((EnergyMax2D[:,2], EnergyMax2D[:,1], EnergyMax2D[:,0]))]
  WpiMax2Dsorted     = WpiMax2D[np.lexsort((WpiMax2D[:,2], WpiMax2D[:,1], WpiMax2D[:,0]))]
  
  del TotalField2D, DensityMax2D, CurrentMax2D, EnergyMax2D, WpiMax2D, DensityMax2D_grad
  
  ## Lazy to build up a new function for that for now.
  unit=1E-9
  PlotGradient=False
  if(WhatToPlot == BicolorE1E2plotDictionnary['Nexc']): 
      QuantityTitle=r'$N_{exc}$ (m$^{-3}$), $(\lambda_1, \lambda_2) = $('+str(round(wavelength1*1E9,0))+', '+str(round(wavelength2*1E9,0))+') nm'
      ExportToTxt(DensityMax2Dsorted, filename+".csv")
  elif(WhatToPlot == BicolorE1E2plotDictionnary['NexcGrad']): #We'll plot grad(Nexc)_Intensity to find number of required photons
      filename=filename_gradNexc
      DensityMax2Dsorted = DensityMax2DgradSorted
      QuantityTitle=r'$N_{exc}$, $(\lambda_1, \lambda_2) = $('+str(round(wavelength1*1E9,0))+', '+str(round(wavelength2*1E9,0))+') nm'
      ExportToTxt(DensityMax2DgradSorted, filename+".csv")
      PlotGradient=True
      unit=1E0 #As we compute gradient, let's keep SI units here. 
  elif(WhatToPlot == BicolorE1E2plotDictionnary['E2']): 
      DensityMax2Dsorted = TotalField2Dsorted; 
      filename = filename_field
      QuantityTitle=r"$E_{tot}^2$ (V/nm)$^2$, $(\lambda_1, \lambda_2) = $("+str(round(wavelength1*1E9,0))+", "+str(round(wavelength2*1E9,0))+") nm."
      ExportToTxt(TotalField2Dsorted, filename+".csv")
  elif(WhatToPlot == BicolorE1E2plotDictionnary['J']): 
      #PlotGradient=False
      DensityMax2Dsorted = CurrentMax2Dsorted; 
      filename = filename_current
      QuantityTitle=r"$J_{end}$ (at.u.), $(\lambda_1, \lambda_2) = $("+str(round(wavelength1*1E9,0))+", "+str(round(wavelength2*1E9,0))+") nm."
      ExportToTxt(CurrentMax2Dsorted, filename+".csv")
  elif(WhatToPlot == BicolorE1E2plotDictionnary['Energy']): 
      #PlotGradient=False
      DensityMax2Dsorted = EnergyMax2Dsorted
      filename = filename_energy
      QuantityTitle=r"$\xi_{end}$ (eV), $(\lambda_1, \lambda_2) = $("+str(round(wavelength1*1E9,0))+", "+str(round(wavelength2*1E9,0))+") nm."
      ExportToTxt(EnergyMax2Dsorted, filename+".csv")
  elif(WhatToPlot == BicolorE1E2plotDictionnary['Wpi']): 
      PlotGradient=True
      DensityMax2Dsorted = WpiMax2Dsorted
      filename = filename_wpi
      QuantityTitle=r'$w_{PI}$ (W$^{-1}$.m$^{-1}$), $(\lambda_1, \lambda_2) = $('+str(round(wavelength1*1E9,0))+', '+str(round(wavelength2*1E9,0))+') nm'
      ExportToTxt(WpiMax2Dsorted, filename+".csv")
  else: 
      PlotGradient=False
      print Header+"** plot_Zhukov_BicolorMaps: wrong request!"
  
  
  # Shaping the tables
  BrutX = unit*DensityMax2Dsorted[:,0] #Efield1
  BrutY = unit*DensityMax2Dsorted[:,1] #Efield2
  BrutZ = DensityMax2Dsorted[:,2]      #Quantity to plot
  
  ### One need to extract interpolation of BrutZ on a regular mapping. 
  if(not len(BrutZ)==0): 
      RegX, RegY, RegZ = plot2dHeatMap(BrutX, BrutY, BrutZ, QuantityTitle, filename, ShowFinalPlot, PlotGradient)
  else: 
      print(Header+"** Error: Selected database is empty.")
      exit()
  
  ## SPECIAL CASE: preparing figure to support Nadya's idea that BICOLOR induces higher absorption than separed single-colors. 
  ## AbsorbedEnergy(Field or intensity)
  ## Intensity1 
  ## Intensity2
  ## Intensity_1+2 ( t )  = 0.5*c*eps0*|E1(t)+E2(t)|**2 #where is optical index? 
  ## We shall select the peak intensity. 
  ## 
  
  if(WhatToPlot == BicolorE1E2plotDictionnary['Energy']):
    # Back to SI unit, but on a regular mesh. 
    Field1   = RegX/unit #back to V/m. 
    Field2   = RegY/unit #back to V/m. 
    
    # 
    
    ## WARNING: induced fields should be necessary for this computation. 
    # Selecting one single line to plot in these data. 
    plt.figure()
    plt.plot()
  
  #if(WhatToPlot == BicolorE1E2plotDictionnary['Nexc']):
    #print "== Performing a section in the two-dimensional Nexc profile."
    ## We will use interp2d in the Brut data
    
    #DB_Field1 = BrutX
    #DB_Field2 = BrutY
    #DB_Wpi    = BrutZ
    
    #print Header+"Data well imported from DB."
    #print DB_Field1.shape, DB_Field2.shape, DB_Wpi.shape
    ### Attempting to use lists of points for interp2d
    
    #WPI_func = interp2d(DB_Field1, DB_Field1, DB_Wpi, fill_value=-1) #interfaces w_PI with E1, E2 values, as a function in
    
    ##EfieldSweep = np.linspace(1E9,1E10, 10)
    #print WPI_func(1E8, 1E8) ##Weirdly returns something out of bounds
    #exit()
    ## Do we need a 2D matrix interpolation? 
    ## This works for 2D matrix interpolation
    #print Header+"** Preparing interpolation of w_pi(E1,E2)..."
    #xdim = int(np.sqrt(len(DB_Field1))) #NOTE: assumes data are square. True for 800x1600, CEP=0 at least.
    #print xdim
    #Wpi_2D_X  = DB_Field1.reshape(xdim, xdim)
    ##print Wpi_2D_X
    #Wpi_2D_Y  = DB_Field2.reshape(xdim, xdim)
    ##print Wpi_2D_Y
    #Wpi_2D    = DB_Wpi.reshape((xdim, xdim))
    ##print Wpi_2D
    #Wpi_X     = Wpi_2D_X[:,0]
    ##print Wpi_X
    #Wpi_Y     = Wpi_2D_Y[0,:]
    ##print Wpi_Y
    #Wpi_2D_t = np.transpose(Wpi_2D)
    ##print Header+"** Interpolating the w_PI..."
    ##InterpolationOrder=1
    ## WARNING: the 800x1600 files have #1: E2**1 !! and #2: E1**2 ! Therefore we had to transpose the matrix. 
    #WPI_func = interp2d(Wpi_X, Wpi_Y, Wpi_2D_t) #interfaces w_PI with E1, E2 values, as a function in RxR. 
    ##WPI_func = interp2d(Wpi_X, Wpi_Y, Wpi_2D)
    #print Header+"Trying to 2D-interpolate at one (E1,E2) point: "
    #print WPI_func(0., 0.) #**2, FieldEnvelopeNormalized2_CGS**2)
  
  del DensityMax2Dsorted
  
  #print "=== CONSTRUCTING the Zhukov maps from the bicolor Keldysh he provided."
  
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
      runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  W_pi, PulseType, RiseTime, FinalTime = ExtractData(PlottingDoublePulses_SelectedColor2)
      
      print "Color1: "+str(Color1)+" m, Color2: "+str(Color2)+" m."
      if(len(Ne_final) >= 1):
        ToBeAdded        = np.array([wavelengths1[0], wavelengths2[0], Ne_final[0]])
        ToBeAddedCurrent = np.array([wavelengths1[0], wavelengths2[0], TotalCurrent_final[0]])
        ToBeAddedEnergy  = np.array([wavelengths1[0], wavelengths2[0], ElectronEnergy_max[0]])
        ToBeAddedWpi     = np.array([wavelengths1[0], wavelengths2[0], W_pi[0]])
        print ToBeAdded, ToBeAddedCurrent, ToBeAddedEnergy, ToBeAddedWpi
        DensityMax2D=np.vstack((DensityMax2D, ToBeAdded))
        CurrentMax2D=np.vstack((CurrentMax2D, ToBeAddedCurrent))
        EnergyMax2D =np.vstack((EnergyMax2D, ToBeAddedEnergy))
        Wpi2D       =np.vstack((Wpi2D, ToBeAddedWpi))
        del ToBeAdded, ToBeAddedCurrent, ToBeAddedEnergy, ToBeAddedWpi
        
  # Stack is built. Export, show and plot on a heat map. 
  print "** Info: exporting the sub-plot2D: Ne_final(wavelength, Field), TotalCurrent_Final(), ElectronEnergy_max()..."
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

def plotAdjustedCollisionFrequencies_wavelength(): #{{{
# Database was constructed by hand

    Efield = np.array(          [1E9,      2E9,       3E9,       4E9,       5E9,      1E10])
    CollTime_1200nm = np.array( [126E-15,  22.75E-15, 16.05e-15, 14.1e-15,  1.27e-14, 1.24e-14])
    CollTime_800nm = np.array(  [4.85e-15, 5.6E-15,   5.9E-15,   6.55E-15,  7.25e-15, 6.9e-15])
    CollTime_600nm = np.array(  [np.nan,   7.3E-15,   6.95E-15,  7.2E-15,   6.85e-15, 5.60e-15])
    CollTime_400nm = np.array(  [1.45e-15, 1.55e-15,  1.75E-15,  2.1E-15,   2.65e-15, 3.55e-15])

    # Select the runs where Collision Frequency is provided. 
    # Plot as function of the peak field amplitude. 

    plt.figure()
    plt.xlabel(r'Field peak amplitude $E_{\mathrm{peak}}$ (V/nm)')
    plt.ylabel(r'Estimated collision time $\nu^{-1}$ (fs)')
    plt.loglog(Efield*1E-9, 1E15*CollTime_400nm,  'b-+',  label=r"$\lambda=400$ nm")
    plt.loglog(Efield*1E-9, 1E15*CollTime_600nm,  'g-+',  label=r"$\lambda=600$ nm")
    plt.loglog(Efield*1E-9, 1E15*CollTime_800nm,  'k-+',  label=r"$\lambda=800$ nm")
    plt.loglog(Efield*1E-9, 1E15*CollTime_1200nm, 'r-+', label=r"$\lambda=1200$ nm") 
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig("CollisionFrequencies.eps")
    #plt.show()
    
#}}}                                                           

## Prepare all plots for single-color excitation with VP Zhukov [Phys. Rev. B]
def plot_Zhukov_SingleColor_LDA(): #{{{
    Egap_SI = 2.56*e
    #Egap_SI = 1.55*e #Dressed in presence of driving photon field
    #Egap_SI = 1.5*e #Manually adjusted gap. 
    Egap_SI_indirect = 0.51*e
    PulseType  = 3 #STH: 2, StH: 3
    wavelength = 800e-9
    ShowPlot   = 0
    PlotWpiAxis = True
    MoveIon=0
    
    ## These are experimental data from Palik. More correct would be to use linear response computed from TDDFT. 
    OpticalIndex={ #
        483e-9: 4.4028, 484e-9: 4.3964,
        800e-9: 3.6924167231905, 1600e-9: 3.4826, 
        2200e-9: 3.4548, 2600e-9: 3.4457, 3200e-9: 3.43438589865999}
    
    ## NOTE: DUMMY OpticalIndex account
    #OpticalIndex={
        #483e-9: 1., 484e-9: 1.,
        #800e-9: 1., 1600e-9: 1., 
        #2200e-9: 1., 2600e-9: 1., 3200e-9: 1. }
    
    #Egap_SI_dressed = {483e-9:, 484e-9: ,
                       #800e-9: 1.55*e, 1600e-9: ,
                       #2200e-9:, 2600e-9: }
    
    ## Plots effect of the STH Rise Time. We even found an incomplete computations thanks to this. Works super well, then. 
    RiseCycles = 0 #0: undefined, mixes everything. Basically, we must redefine all runs which had no well defined RiseCycle so far. 
    Efield     = 1E9
    #plotEffectOfRiseTime(PlottingTable, wavelength, 20e-15, Efield, PulseType, ShowPlot, Egap_SI)
    
    # This is a manually prepared figure. 

    plotAdjustedCollisionFrequencies_wavelength()
    
    ## Effects of CEP on single pulse excitation, current and energy, and harmonic spectrum
    RiseCycles = 2
    Efield     = 1E9
    # plotEffectOfCEP(PlottingTable, wavelength, 20e-15, Efield, PulseType, RiseCycles, ShowPlot, Egap_SI)
    #plotEffectOfCEP(PlottingTable, wavelength, 10e-15, Efield, PulseType, RiseCycles, ShowPlot)
    Efield     = 5E9 #NOTE: Saturation may prevent to see any difference. 
    # plotEffectOfCEP(PlottingTable, wavelength, 20e-15, Efield, PulseType, RiseCycles, ShowPlot, Egap_SI)
    #plotEffectOfCEP(PlottingTable, wavelength, 10e-15, Efield, PulseType, RiseCycles, ShowPlot)
    #NOTE: === LOCAL FILTERING
    #PlottingTable_filtered = FilterDatabaseRemove(PlottingTable, 30e-15, IndexTable['pulses1'])
    PlottingTable_filtered = PlottingTable
    

    ## Effect of intensity figures - not much interesting for "long" pulses
    # plotEffectOfIntensity_CEP(PlottingTable_filtered, wavelength, 20e-15, False, 0.0, PulseType, RiseCycles, ShowPlot, Egap_SI)
    #plotEffectOfIntensity_CEP(PlottingTable_filtered, wavelength, 10e-15, False, 0.0, PulseType, RiseCycles, ShowPlot)
    #plotEffectOfIntensity_CEP(PlottingTable_filtered, wavelength, 30e-15, False, 0.0, PulseType, RiseCycles, ShowPlot)
    
    ## ABSORPTION PROBABILITY FOR 800 nm - Phys Rev B Derrien et al 2020
    plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI, Egap_SI_indirect, OpticalIndex[wavelength]) #This one construct the ionization probabilities! 
    
    ## RESONANT ABSORPTION - not published, not understood
    wavelength = 483e-9 #Resonant wavelength results for LDA.
    #plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI, Egap_SI_indirect, OpticalIndex[wavelength]) #This one construct the ionization probabilities! 
    wavelength = 484e-9 #Resonant wavelength results for LDA.
    #plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI, Egap_SI_indirect, OpticalIndex[wavelength]) #This one construct the ionization probabilities! 
    
    ## READY
    wavelength = 800e-9
    #plotEffectOfIntensity_wavelengths(PlottingTable_filtered, 20e-15, 0., False, wavelength, PulseType, RiseCycles, ShowPlot, Egap_SI) 
    ##plotEffectOfIntensity_wavelengths(PlottingTable_filtered, 10e-15, 0., False, wavelength, PulseType, RiseCycles, ShowPlot)
    #plotEffectOfIntensity_wavelengths(PlottingTable_filtered, 30e-15, 0., False, wavelength, PulseType, RiseCycles, ShowPlot, Egap_SI)
    
    ## CHECKING DRUDE MODEL VALIDITY (results as function of wavelength)
    ShowPlot=1
    #plotEffectOfWavelength_2dMap(PlottingTable, 20e-15, PulseType, RiseCycles, ShowPlot)
    plotEffectOfWavelength_Tolerance(PlottingTable, 1E-4, 20e-15, 0., PulseType, RiseCycles, ShowPlot, Egap_SI, Egap_SI_indirect, 500e-9, 3300e-9)
    
    ## Checking what is happening in the mid-IR world. NOTE: NOT ENOUGH DATA? 
    #plotEffectOfWavelength_Tolerance(PlottingTable, 1E-4, 60e-15, 0., PulseType, RiseCycles, ShowPlot, Egap_SI, Egap_SI_indirect, 3e-6, 25e-6)
    #plotEffectOfWavelength_Tolerance(PlottingTable, 0.05, 30e-15, 0., PulseType, RiseCycles, ShowPlot, Egap_SI, Egap_SI_indirect) #ALMOST NO DATA for this one

    ### OLD STUFFS
    #plotEffectOfIntensity_CEP(PlottingTable, wavelength, 10e-15, PulseType, RiseCycles, ShowPlot)
    ## Plots the N_exc(wavelength, field) (2D plot) # NOTE: this can be compared with GrojoPRB2013. 
    
    ## MID-IR: Preparation of the mid-IR experiment: measuring the cross section absorption in mid-IR 
    #wavelength = 1600e-9
    #plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, 1600e-9, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI, Egap_SI_indirect, OpticalIndex[wavelength]) #contains MoveIon, yes?
    
    # Some data are available for 2200 nm, but they are very sparse (not priority)
    #plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, 2200e-9, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI, Egap_SI_indirect, OpticalIndex[wavelength]) #contains MoveIon, yes? 
    
    #plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, 2600e-9, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI, Egap_SI_indirect, OpticalIndex[wavelength]) #contains MoveIon, yes?
    #plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, 3200e-9, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI, Egap_SI_indirect) #contains MoveIon, yes? 
    
    ## We now control in presence in MoveIons
    try:
        wavelength = 800e-9
        plotEffectOfIntensity_PulseDuration(PlottingTable_filtered,  wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI, Egap_SI_indirect, OpticalIndex[wavelength], MoveIon)
    except:
        print Header+"** Warning: no data could be found, or bug took place for 800 nm."
    try:
        wavelength = 1600e-9
        plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI, Egap_SI_indirect, OpticalIndex[wavelength], MoveIon)
    except:
        print Header+"** Warning: no data could be found, or bug took place for 1600 nm."
    #try:
    wavelength = 3200e-9
    plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI, Egap_SI_indirect, OpticalIndex[wavelength], MoveIon)
    #except:
        #print Header+"** Warning: no data could be found, or bug took place for 3200 nm."
#}}}
def plot_Zhukov_SingleColor_TB09(): #{{{
    Egap_SI = 3.04*e
    Egap_SI_indirect = 0.98*e
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
    #ShowPlot=1
    plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI, Egap_SI_indirect)
    
    #This one construct the ionization probabilities! 
    #wavelength = 483e-9 #Resonant wavelength results for LDA.
    wavelength = 407.58e-9 #Resonant wavelength results for TB09. 
    plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI, Egap_SI_indirect) #This one construct the ionization probabilities! 
    #exit()
    #wavelength = 484e-9 #Resonant wavelength results for LDA.
    #plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot) #This one construct the ionization probabilities! 
    ## READY
    wavelength = 800e-9
    plotEffectOfIntensity_wavelengths(PlottingTable_filtered, 20e-15, 0., PlotWpiAxis, wavelength, PulseType, RiseCycles, ShowPlot, Egap_SI) 
    #plotEffectOfIntensity_wavelengths(PlottingTable_filtered, 10e-15, 0., PlotWpiAxis, wavelength, PulseType, RiseCycles, ShowPlot)
    #plotEffectOfIntensity_wavelengths(PlottingTable_filtered, 30e-15, 0., PlotWpiAxis, wavelength, PulseType, RiseCycles, ShowPlot)
    
    ShowPlot=1
    #plotEffectOfWavelength_2dMap(PlottingTable, 20e-15, PulseType, RiseCycles, ShowPlot)
    plotEffectOfWavelength_Tolerance(PlottingTable, 0.03, 20e-15, 0., PulseType, RiseCycles, ShowPlot, Egap_SI, Egap_SI_indirect)
    #plotEffectOfWavelength_Tolerance(PlottingTable, 0.03, 30e-15, 0., PulseType, RiseCycles, ShowPlot, Egap_SI)
    ### OLD STUFFS
    #plotEffectOfIntensity_CEP(PlottingTable, wavelength, 10e-15, PulseType, RiseCycles, ShowPlot)
    ## Plots the N_exc(wavelength, field) (2D plot) #NOTE: this can be compared with GrojoPRB2013. 

#}}}

## Prepare all plots for bicolor paper with V. Zhukov [1 Phys. Rev. Lett + 1 Phys. Rev. A]
def plot_Zhukov_BicolorMaps(): #{{{
    Egap_SI = 2.56*e; Egap_SI_indirect=0.51*e
    
    ## Plotting E1, E2 maps 
    for tau in [30e-15]:
        PulseType  = 3; RiseCycles = 2 #StH corrected. 
        #PulseType  = 2;  RiseCycles = 0 #bugged STH. 
        
        ##Bicolor_plotEffectOfCEP(PlottingTable, 1030e-9, tau, 1E9, 515e-9, tau, 1E9, PulseType, RiseCycles)
        ##Bicolor_plotEffectOfCEP(PlottingTable, 1030e-9, tau, 9E9, 515e-9, tau, 9E9, PulseType, RiseCycles)
        ##exit()
        
        ## DONE: 800+1600
        ## TODO: 1030+515, 800+400
        ## TODO: 800+1030
        
        #### These runs are separatly SUB-band gap, but jointly ABOVE band gap. 
        ### Commensurable 800+1600
        #wavelength1 = 800e-9; wavelength2 = 1600e-9
        CEP1=0.; CEP2=0.; Egap=2.56*e; meff=0.2226; 
        ## BicolorE1E2plotDictionnary={'Nexc': 0, 'E2': 4, 'J': 1, 'Energy': 2, 'Wpi': 3}
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['Nexc'])
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['NexcGrad'])
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['J'])
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['Energy'])
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['Wpi'])
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['E2'])
        ### Generate the tables from Zhukov files and plot them
        ##VPZ_Wpi2D(wavelength1, wavelength2, 4E9, 2E9)
        
        ##VZ_generateWpiTables(Efield1, Efield2, wavelength1, wavelength2, CEP1, CEP2, Egap, meff)
        ##plot_VZ_Bicolor()
        
        ### Commensurable ## NOTE: these results will never reveal multiphotonic steps, because of free-carrier absorption, but also because of the multiple bands. Only possibility is to make sections for a given field, and to analyze the number of multiphoton absorption. 
        ## BUG: some new runs on Prometheus keep RiseCycle=0 although they should be RiseCycle=2. 
        ## Then, I impose RiseCycles = 0 to avoid selecting them. TODO: fix it when duplicates with various RISE cycles would appear. 
        wavelength1 = 1030e-9; wavelength2 = 515e-9; RiseCycles=0 #NOTE: just to lift selection
        plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['Nexc'])
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['NexcGrad'])
        plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['J'])
        plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['Energy'])
        plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['Wpi'])
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['E2'])
        ##VZ_generateWpiTables(Efield1, Efield2, wavelength1, wavelength2, CEP1, CEP2, Egap, meff)
        ##VPZ_Wpi2D(wavelength1, wavelength2, 0.8E10, 1E10)
        ### Incommensurable
        #wavelength1 = 800e-9 
        #wavelength2 = 1030e-9; 
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  0)
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  1)
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  2)
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  3)
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  4)
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  5)
                
        ### Plotting the Stark effect for 2 colors (1: produces N_exc(omega). 2: Energy(omega).)
        
    # Plotting f(wavelength2) for several fields, and a given wavelength1    
    for tau in [30e-15]: #{{{
        ShowPlot=1
        FieldLimit=1E10 #5E9
        wavelength1 = 1030e-9;
        PulseType  = 3; RiseCycles = 0 #NOTE: RiseCycles=2 would be better #StH corrected. 
        for Field1 in [2.5E8]: #, 1.0E9]: # TODO: how to add the E2=0 & E1!=0 ?
            #try:
            Bicolor_plotEffectOfWavelength_Tolerance(PlottingTable, wavelength1, 1E-5, Field1, FieldLimit, tau, 0., PulseType, RiseCycles, ShowPlot, Egap_SI, Egap_SI_indirect)
            print(Header+"Bicolor_plotEffectOfWavelength_Tolerance() success for (tau, wavelength, Efield1)=("+str(tau*1E15)+" fs, "+str(wavelength1*1E9)+" nm, "+str(Field1*1E-9)+" V/nm)")
            #except:
                #print(Header+"Bicolor_plotEffectOfWavelength_Tolerance() failed for (tau, wavelength, Efield1)=("+str(tau*1E15)+" fs, "+str(wavelength1*1E9)+" nm, "+str(Field1*1E9)+" V/nm)")
    exit()
    #}}}
    #VPZ_Wpi2D(wavelength1, wavelength2, 0.8E10, 1E10)
    
    ## Stark analysis of the bicolor simulations. 
    # Plot N_exc(E2; E1=fixed)
    
    
    # How to plot absorbed energy as function of the CEP2. 
#}}}

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
    
def plot_HHG():
        for tau in [20e-15, 30e-15, 50e-15]:
            PulseType = 3 #1: SinSquared. 3: STH
            RiseCycles = 0 #We dont care
            wavelength1 = 3100e-9
            ## Which type of plot for HHG? 
            # HHG(omega; run1, run2, run3, run4, run5) #to check convergency [USEFUL]
            # HHG(omega; LAMBDA1, LAMBDA2, given fields) [RISKY, low visibility]
            # HHG(omega; E1, E2; LAMBDA1, LAMBDA2) [EFFICIENT, VISIBLE, MEANINGFUL] Expected: increase of the piedestal/cutoff
            # HHG(omega; tau1; LAMBDA1) Expected: reduction of the harmonics with time. Gabor transform necessary. 
            # HHG(omega; lambda1) #Figure made by N. Tancogne in his PRL 2016. Great to see the cutoff. 
            # Filter: we don't care about CEP
            
    

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
#plot_Zhukov_SingleColor_TB09()

## Zhukov bicolor plotting
plot_Zhukov_BicolorMaps() #LDA gap

## Rubio shift current plotting
#plot_ShiftCurrents()

## Preparing HHG spectra dependencies for Hamed Merdji
#plot_HHG()

print Header+"End of program."
exit()

