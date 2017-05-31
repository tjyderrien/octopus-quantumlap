#!/usr/bin/env python
#-*- coding: utf-8 -*-

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
import glob

import calculateVolume as Sample
# import libKeldysh as Keldysh
from libAtomicUnits import *
from libPlottingSettings import *
from libString import *
from libPlotting import *
from libKeldysh import *
import libDatabase

Header="[plotFinalQuantities] "

ShowFinalPlot              = 1

EnableDoublePulseAnalysis  = 2 #1: yes, 2: auto on command line
EnableInducedFieldAnalysis = 2 #1: yes, 2: auto on command line

Density=5E28

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
path='*/*/output_iter/' #valid for the runs computed from a same k-grid folder. 
DatabaseName = 'Database.csv'

print Header+"Attempting to load previous results..."
try:
  PlottingTable = loadtxt(DatabaseName, delimiter="\t")
except:
  PlottingTable = ListAllOutput(path)

# Save the table to a file
libDatabase.ExportToTxt(PlottingTable, DatabaseName)

def ExtractData(PlottingTable): #{{{
  wavelengths1      = PlottingTable[:,0]
  pulses1           = PlottingTable[:,1]
  Efield1           = PlottingTable[:,2]
  CEP1              = PlottingTable[:,3]

  wavelengths2      = PlottingTable[:,4]
  pulses2           = PlottingTable[:,5]
  Efield2           = PlottingTable[:,6]
  CEP2              = PlottingTable[:,7]

  Delay             = PlottingTable[:,8]
  CrossPolarization = PlottingTable[:,9]

  #TheoryLevel      = PlottingTable[:,10]
  FieldInduced      = PlottingTable[:,10]

  Ne_max            = PlottingTable[:,11]
  Ne_final          = PlottingTable[:,12]
  irreversibility   = PlottingTable[:,13]
  Energy_max        = PlottingTable[:,14]

  return wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, FieldInduced, Ne_max, Ne_final, irreversibility, Energy_max
#}}}

wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, FieldInduced, Ne_max, Ne_final, irreversibility, Energy_max = ExtractData(PlottingTable)
#TODO: add this field ! #TheoryLevel,

IndexTable = {'wavelengths1': 0, 'pulses1': 1, 'Efield1': 2, 'CEP1': 3, 'wavelengths2': 4, 'pulses2': 5, 'Efield2': 6, 'CEP2': 7, 'Delay':  8, 'CrossPolarization': 9, 'FieldInduced': 10, 'Ne_max': 11, 'Ne_final': 12, 'irreversibility': 13, 'Energy_max': 14}
#TODO: add 'TheoryLevel': 10, 

print Header+"** Preparing the summary plots..."

def CountTheNumberOfLaserParameters(Database):#{{{
  wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, FieldInduced, Ne_max, Ne_final, irreversibility, Energy_max = ExtractData(Database)
  # Generates lists of unique entries to count the range of parameters we have. 
  wavelengths1_uniq      = list(set(wavelengths1))
  wavelengths2_uniq      = list(set(wavelengths2))
  Efield1_uniq           = list(set(Efield1))
  Efield2_uniq           = list(set(Efield2))
  pulses1_uniq           = list(set(pulses1))
  pulses2_uniq           = list(set(pulses2))
  Delay_uniq             = list(set(Delay))
  CrossPolarization_uniq = list(set(CrossPolarization))
  CEP1_uniq              = list(set(CEP1))
  CEP2_uniq              = list(set(CEP2))

  wavelength_numbers     = len(wavelengths1_uniq)
  wavelength2_numbers    = len(wavelengths2_uniq)


  print "====== PULSE 1 ====="
  print "** Efield1*           size = "+str(len(Efield1_uniq))
  print "** wavelength1*       size = "+str(len(wavelengths1_uniq))
  print "** pulse duration*    size = "+str(len(pulses1_uniq))
  print "** CEP1               size = "+str(len(CEP1_uniq))

  print "====== PULSE 2 ====="
  print "** Efield2*           size = "+str(len(Efield2_uniq))
  print "** wavelength2*       size = "+str(len(wavelengths2_uniq))
  print "** pulse duration*    size = "+str(len(pulses2_uniq))
  print "** CEP2               size = "+str(len(CEP2_uniq))

  print "=== COUPLING EFFECTS ==="
  print "** Delays             size = "+str(len(Delay_uniq))
  print "** CrossPolarization  size = "+str(len(CrossPolarization_uniq))

  print "====== RESULTS ====="
  print "** Ne_final           size = "+str(len(Ne_final))
  print "** Irreversibility    size = "+str(len(irreversibility))
  print "** Total energy       size = "+str(len(Energy_max))
  
  return(wavelengths1_uniq, CEP1_uniq)
#}}}

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
SizeX = 10
SizeY = 10

plt.figure(figsize=(SizeX,SizeY))

## We change strategy: now we plot all the available data based on the number of parameters we have detected. The input parameters will not define what we have to plot anymore. We prepare all the plots, which will simplify the understanding of this file. 

#====== Plot: SINGLE PULSE, SAME WAVELENGTHS, SAME POLARIZATIONS: final density as function of field ======
#Filter: double pulses, same wavelengths, same polarization
PlottingSinglePulses = PlottingTable #get the full database
PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0. 
PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['CEP1']) #removes fields where CEP != 0. 
#PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 800e-9, IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts.

## prints the current state of selection
wavelength_uniq, CEP1_uniq = CountTheNumberOfLaserParameters(PlottingSinglePulses)

# Simple pulse, single color

if(len(wavelength_uniq) == 1): #{{{ #single color simulations
  # Complexity level: 1
  
  ax1 = plt.subplot(plotA)
  
  pplt.title('Density as function of peak intensity')
  # Number of sweeped parameters: 3 (Efield, Delay, CrossPolarization)
  #pplt.scatter(FieldToIntensity(Efield1), Ne_final, c=Delay, s=(30.+CrossPolarization*50.), label=r'$N_{exc}$[Field]')
  for wavelength_value in wavelength_uniq:
    # select the corresponding wavelength in results
    PlottingSinglePulses_wavelength = FilterDatabase(PlottingSinglePulses, 800e-9, IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts.

    wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, FieldInduced, Ne_max, Ne_final, irreversibility, Energy_max = ExtractData(PlottingSinglePulses_wavelength)

    pplt.plot(FieldToIntensity(Efield1), Ne_final, 'o', label=r'$N_{exc}$, $\lambda=$'+str(round(wavelength_value*1E9))+' nm')
  ax1.legend(loc=2, borderaxespad=0.) #bbox_to_anchor=(1.05, 1), 
  #pplt.xlabel(r'Peak intensity (W/m$^2$)')
  pplt.ylabel(r'Final density (m$^{-3}$)')
#  pplt.ylim((1E16,4.*5e28))
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  
  ax2 = plt.subplot(plotB)
  plt.title('Irreversibility factor')
  pplt.plot(FieldToIntensity(Efield1), irreversibility, 'o', label='Irreversibility')
  pplt.xlabel(r'Intensity (W/m$^2$)')
  pplt.ylabel('Irreversibility factor')
  pplt.ylim((0,1))
  ax2.set_xscale('log')
  ax22 = ax2.twinx()
  for wavelength_value in wavelength_uniq:
    # select the corresponding wavelength in results
    PlottingSinglePulses_wavelength = FilterDatabase(PlottingSinglePulses, 800e-9, IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts.

    wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, FieldInduced, Ne_max, Ne_final, irreversibility, Energy_max = ExtractData(PlottingSinglePulses_wavelength)
    pplt.plot(FieldToIntensity(Efield1), gammaKeldysh(2.56*e, 0.228, Efield1, wavelengths1), 'or', label=r'Adiabaticity, $\lambda=$'+str(round(wavelength_value*1E9))+' nm')
  pplt.ylabel('Adiabadicity factor')
  #pplt.ylim((0,1))
  ax22.set_yscale('log')
  
  # ax2.set_yscale("log")
  ax2.legend(loc=2, borderaxespad=0.)
  # ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
  ax22.legend(loc=1, borderaxespad=0.)
  
  # Saving the plot
  pplt.savefig('Metaplot-SinglePulse-SingleColor-EffectOfWavelength.eps')
  pplt.savefig('Metaplot-SinglePulse-SingleColor-EffectOfWavelength.png')
  #if(ShowFinalPlot==1):
    #pplt.show()

del ax1, ax2, ax22, PlottingSinglePulses, wavelength_uniq, CEP1_uniq
#}}}
  
#====== Plot: SINGLE PULSE, EFFET OF THE CEP, SAME WAVELENGTHS, SAME POLARIZATIONS: final density as function of field ======
#Filter: double pulses, same wavelengths, same polarization

# Simple pulse, single color
PlottingSinglePulses = PlottingTable #get the full database
PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0. 
PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 800e-9, IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts. #TODO: this field should be making a loop in available wavelengths given in wavelengths1_uniq. 

wavelength_uniq, CEP1_uniq = CountTheNumberOfLaserParameters(PlottingSinglePulses)

#if(len(CEP1_uniq) == 1): #{{{ #single color simulations
# Complexity level: 2

ax1 = plt.subplot(111)

pplt.title('Final density as function of peak intensity, for various CEP')
# Number of sweeped parameters: 3 (Efield, Delay, CrossPolarization)
#pplt.scatter(FieldToIntensity(Efield1), Ne_final, c=Delay, s=(30.+CrossPolarization*50.), label=r'$N_{exc}$[Field]')
for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
  PlottingSinglePulses_CEP = FilterDatabase(PlottingSinglePulses, CEP1_value, IndexTable['CEP1']) #select the local CEP
  ## prints the current state of selection and stores the values of sweeped parameters
  wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, FieldInduced, Ne_max, Ne_final, irreversibility, Energy_max = ExtractData(PlottingSinglePulses_CEP)
  
  pplt.plot(FieldToIntensity(Efield1), Ne_final, 'o', label=r'$N_{exc}^{final}$, CEP='+str(round(CEP1_value*180./pi)))
    
ax1.legend(loc=4, borderaxespad=0.) #bbox_to_anchor=(1.05, 1), 
pplt.xlabel(r'Peak intensity (W/m$^2$)')
pplt.ylabel(r'Final density (m$^{-3}$)')
#  pplt.ylim((1E16,4.*5e28))
ax1.set_xscale('log')
ax1.set_yscale('log')

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
pplt.savefig('Metaplot-SinglePulse-SingleColor-EffectOfCEP.eps')
pplt.savefig('Metaplot-SinglePulse-SingleColor-EffectOfCEP.png')
#if(ShowFinalPlot==1):
  #pplt.show()

del ax1, PlottingSinglePulses
#}}}

exit()

##====== Plot: DOUBLE PULSE, SAME WAVELENGTH, SAME POLARIZATION: final density as function of delay and field ======
  
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
  
 

##TODO: shall we make a multi-wavelength figure?
##TODO: - If wavelength2!=wavelength1, yes.
##TODO: - If len(set(wavelength1)) > 1, yes.  

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

#====== Plot: DOUBLE PULSE, SAME WAVELENGTH, SAME POLARIZATION: final density as function of delay and field ======
#Filter: double pulses, same wavelengths, same polarization
PlottingDoublePulses = PlottingTable #get the full database
PlottingDoublePulses = FilterDatabaseGreaterThan(PlottingDoublePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak < 1 V/m. #works
PlottingDoublePulses = FilterDatabase(PlottingDoublePulses, wavelengths1, IndexTable['wavelengths2']) #removes fields where wavelength2 != wavelength1
PlottingDoublePulses = FilterDatabase(PlottingDoublePulses, float(0.), IndexTable['CrossPolarization']) #removes fields where wavelength2 != wavelength1
wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, FieldInduced, Ne_max, Ne_final, irreversibility, Energy_max = ExtractData(PlottingDoublePulses)

#ax1 = plt.subplot(plotA)

# Number of sweeped parameters: 2
plt.title('Final density (peak intensity, delay)')
print "** 3d plot: Delay*       size = "+str(len(Delay_uniq))
print "** 3d plot: wavelength1* size = "+str(len(wavelengths1_uniq))
#  Ne_final_3d=np.reshape(Ne_final, (len(Efield1_uniq), len(Delay_uniq)))
# Ne_final_3d_interp = interpolate.interp2d(Efield1, Delay, Ne_final)
# Efield1_interp = np.arange(1E8, 1E10, 1E8)
# Delay_interp   = np.arange(0e0, 50e-15, 10e-15)
# Ne_final_3d = Ne_final_3d_interp(Efield1_interp, Delay_interp) 
# plt.contourf(Efield1_interp, Delay_interp, Ne_final_3d, tri=True)
pplt.ylabel('Delay (fs)')
pplt.xlabel('Field ampl. (V/m)')
#pplt.scatter(FieldToIntensity(Efield1), Delay*1E15, s=34+30*Ne_final/Density, c=Ne_final/Density)
pplt.scatter(FieldToIntensity(Efield1), Delay*1E15, c=Ne_final/Density)
pplt.savefig('Metaplot-DoublePulse.eps')
pplt.savefig('Metaplot-DoublePulse.png')
pplt.show()

exit()

