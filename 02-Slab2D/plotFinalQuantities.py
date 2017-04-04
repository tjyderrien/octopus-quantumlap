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

PlottingTable = np.empty((0,7)) 

ShowFinalPlot              = 1

EnableDoublePulseAnalysis  = 2 #1: yes, 2: auto on command line
EnableInducedFieldAnalysis = 2 #1: yes, 2: auto on command line

Density=5E28

if(EnableDoublePulseAnalysis==2 or EnableInducedFieldAnalysis==2): #automatic mode
  try: 
    EnableDoublePulseAnalysis = int(float(sys.argv[1]))
    EnableInducedFieldAnalysis = int(float(sys.argv[2]))
    print "Double pulse mode?"+str(EnableDoublePulseAnalysis)
    print "Induced fields   ?"+str(EnableInducedFieldAnalysis)
#    ShowFinalPlot = 0 #avoids problems with serial plotting 
  except: 
    print "Usage: plotFinalQuantities.py DoublePulse? InducedFields? "
    print "Ex:    plotFinalQuantities.py 1 1"
    exit()
  


# List all output_iter files
for FolderName in glob.glob('*/*/output_iter/'): 
  #TODO: this is not general enough, directory level does not help
  print FolderName

  # Import the quantities from the subfolder
  geometry, Ne, laser, TotalCurrent, Energy = LoadRunOutputs(FolderName)
  # Collect the laser parameters
  filename = FolderName+"../inp"
  Efield2_inp = LookUpStringInFile(filename, "EfieldBsi =")
  Efield2_SI  = float(Efield2_inp[11:15]) #If non-zero, we enable double-pulse analysis. Else: disabled. 
  wavelength_inp = LookUpStringInFile(filename, "lambdaA =")
  wavelength_SI = float(wavelength_inp[10:17])
  PulseDuration_inp = LookUpStringInFile(filename, "tauAsi = ")
  PulseDuration_SI = float(PulseDuration_inp[8:15])
  DoublePulseDelay_inp = LookUpStringInFile(filename, "pulseDelaySI = ")
  DoublePulseDelay_SI = float(DoublePulseDelay_inp[15:22])
  FieldInduced_inp = LookUpStringInFile(filename, "GaugeFieldPropagate = ")
  try:
    FieldInduced_string = FieldInduced_inp[0]
  except:
    FieldInduced_string = ''

  # Filter the induced fields #TODO: reshape these tests. Maybe invert test on InducedFieldDetection and InducedFieldAnalysis?  Problem is: we would like to miss none of the runs. Because counting dots and number of runs is complicated. 
  ## Remove contradictory caes
  if(EnableInducedFieldAnalysis == 1 and FieldInduced_string == "#"):
    print "Jumped: InducedFieldON but DISABLED in computation. [Case 1]"
    continue
  if(EnableInducedFieldAnalysis == 0 and FieldInduced_string == "G" and FieldInduced_inp[22] == "y"):
    print "Jumped InducedFieldOFF but ENABLED in computation. [Case 2]"
    continue
  ## Matching cases also may have some exceptions
  else: 
    if(FieldInduced_string == "#"):
      print "Induced Field was DISABLED in computation."
    elif(FieldInduced_string == "G"):
      print "Induced field was ENABLED in computation."
    else:
      if(EnableInducedFieldAnalysis==1):
        print "Induced field: contradiction. GaugeFieldPropagate was not found in this computation. Jumped."
        continue
      else:
        print "Induced field was not found in this computation. But we take it into account."
  

  # Jump selections which don't match with desire. 
  if(  Efield2_SI > 1E0 and EnableDoublePulseAnalysis == 0): #In case we found a double-pulse
    print "Jump [Case 3]"
    continue #we jump this run, then 
  elif(Efield2_SI < 1E0 and EnableDoublePulseAnalysis == 1): 
    print "Jump [Case 4]"
    continue

#  else: 
#    print "SINGLE-PULSE MODE"
    # Remove the pulses with double

  print "Wavelength (nm): "+str(wavelength_SI)
  print "Pulse duration (s): "+str(PulseDuration_SI)
  print "Double pulse delay (s): "+str(DoublePulseDelay_SI)
  
  # Catch the peak field #We could also take the one of input file..
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
  #print "** Sample volumes **"
  #print "Super-cell total volume = "+str(SuperCellVolume_AU)+" Bohr**3."
  #print "Super-cell total volume = "+str(SuperCellVolume_SI)+" m^{3}." 
  #print ""
  #print "Sample volume = "+str(SampleCellVolume_AU)+" Bohr**3."
  #print "Sample vacuum volume = "+str(SampleCellVolume_SI)+" m^{3}."
  #print ""
  #print "Cell volume = "+str(CellVolume_AU)+" Bohr**3."
  #print "Cell volume = "+str(Volume_AU_to_SI(CellVolume_AU))+" m^{-3}."
  #print "" 
  
  # Catch the final excited electron density. 
  timeNe = np.asfarray(Ne[:,1])
  Nproj_AU = np.asfarray(Ne[:,2]) #number of particles, not a density!
  density_proj_AU = Nproj_AU / SampleCellVolume_AU #This is now a density (a.u.)
  density_proj_SI = Density_AU_to_SI(density_proj_AU)
  sizeOfDensityTab = density_proj_SI.size

  density_proj_SI_final = density_proj_SI[sizeOfDensityTab-1]
  density_proj_SI_max   = density_proj_SI.max()
  density_irreversibility_ratio = density_proj_SI_final / density_proj_SI_max
  
  print "** Electron density analysis "
  print "Final electron density = "+str(density_proj_SI_final)
  print "Max.  electron density = "+str(density_proj_SI_max)
  print "Irreversibility ratio (Final / Max) = "+str(density_irreversibility_ratio)
  ##==== Prepare plotting of energies ====
  timeEnergies   = np.asfarray(Energy[:,1]) #atomic units (! 1 is 2nd column)
  ElectronEnergy = np.asfarray(Energy[:,5]) #Hartree
  TotalEnergy    = np.asfarray(Energy[:,2]) #Hartree
  IonEnergy      = np.asfarray(Energy[:,3]) #Hartree
  Hartree        = np.asfarray(Energy[:,7]) #Hartree
  
  # Append to the table
  # Structure: 
  # wavelength1, wavelength2, tau1, tau2, 
  # Epeak1, Epeak2, 
  # CEP1, CEP2, delay, 
  # kGrid, MeshSize, 
  # NeMax, NeEnd, density_irreversibility_ratio, EnergyMax
  # One line is one octopus run. O_o
  #PlottingTable = np.vstack(( PlottingTable, [ 
				#wavelength1, wavelength2, 
				#tau1, tau2, Epeak1, Epeak2, 
				#CEP1, CEP2, delay, 
				#kGrid, MeshSize, 
				#NeMax, NeEnd, density_irreversibility_ratio, EnergyMax ] ))
  # Simplified version first
  PlottingTable = np.vstack(( PlottingTable, [ 
				wavelength_SI, PulseDuration_SI, Efield_amplitude_SI, 
				DoublePulseDelay_SI, density_proj_SI_max, density_proj_SI_final, 
				density_irreversibility_ratio ] )) ##7
  
  # Cleaning the guys before next loop
  del time, A1x, A2x, timestepTable, Efield, Efield2, geometry, Ne, laser, TotalCurrent, Energy
  #del wavelength1, wavelength2
  #del tau1, tau2
  #del Epeak1, Epeak2
  #del CEP1, CEP2
  #del delay, kGrid, MeshSize
  del density_proj_SI_max, density_proj_SI_final, density_irreversibility_ratio
  #del EnergyMax
  del timeEnergies, ElectronEnergy, TotalEnergy, IonEnergy, Hartree
  

wavelengths1    = PlottingTable[:,0]
pulses1         = PlottingTable[:,1]
Efield1         = PlottingTable[:,2]
Delay           = PlottingTable[:,3]
Ne_max          = PlottingTable[:,4]
Ne_final        = PlottingTable[:,5]
irreversibility = PlottingTable[:,6]

#print wavelengths1.shape, pulses1.shape, Efield1.shape, Delay.shape, Ne_max.shape, Ne_final.shape, irreversibility.shape

## Prepare the summary plots

### Defining plotting structure

#### Plotting: parameters for controlling the canvas
if(EnableDoublePulseAnalysis==0):
  # 2 plots
  plotA=211; plotB=212; #plotC=313
else: #if(EnableDoublePulseAnalysis==1):
# 3 plots
# plotA=311; plotB=312; plotC=313;
  # 4 plots 
  plotA=221; plotB=222; plotC=223; plotD=224;

# Size of main canvas
SizeX = 10
SizeY = 15

plt.figure(figsize=(SizeX,SizeY))

#Plot 1: Final density as function of delay
ax1 = plt.subplot(plotA)
if(EnableDoublePulseAnalysis==1): 
  plt.title('Final density (peak field, delay)')
  print "** 3d plot: Ne_final size = "+str(len(Ne_final))
  Efield1_uniq=list(set(Efield1)); Delay_uniq=list(set(Delay));
  print "** 3d plot: Delay*    size = "+str(len(Delay_uniq)) 
  print "** 3d plot: Efield1*  size = "+str(len(Efield1_uniq))
#  Ne_final_3d=np.reshape(Ne_final, (len(Efield1_uniq), len(Delay_uniq)))
  # Ne_final_3d_interp = interpolate.interp2d(Efield1, Delay, Ne_final)
  # Efield1_interp = np.arange(1E8, 1E10, 1E8)
  # Delay_interp   = np.arange(0e0, 50e-15, 10e-15)
  # Ne_final_3d = Ne_final_3d_interp(Efield1_interp, Delay_interp) 
  # plt.contourf(Efield1_interp, Delay_interp, Ne_final_3d, tri=True)
  pplt.ylabel('Delay (fs)')
  pplt.xlabel('Field ampl. (V/m)')
  pplt.scatter(Efield1, Delay*1E15, s=34+30*Ne_final/Density, c=Ne_final/Density)
else: 
  pplt.title('Density as function of peak field')
  pplt.scatter(Efield1, Ne_final, c=Delay, label=r'$N_{exc}$[Field]')
  # ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
  pplt.xlabel('Field peak amplitude (V/m)')
  pplt.ylabel(r'Final density (m$^{-3}$)')
  pplt.ylim((1E16,4.*5e28))
#  ax1.set_yscale('log')

#Plot 2: irreversibility
ax2 = plt.subplot(plotB)
plt.title('Irreversibility factor')
if(EnableDoublePulseAnalysis==1):
  pplt.ylabel('Delay (fs)')
  pplt.xlabel('Field ampl. (V/m)')
#  area = np.po
  pplt.scatter(Efield1, Delay*1E15, s=34+32*irreversibility, c=irreversibility)
else:
  pplt.scatter(Efield1, irreversibility, c=Delay*1E15, label='Irreversibility')
  # ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
  pplt.xlabel('Field peak amplitude (V/m)')
  pplt.ylabel('Irreversibility factor')
  pplt.ylim((0,1))

if(EnableDoublePulseAnalysis==1):
  #Plot 3: Final density as function of delay
  ax3 = plt.subplot(plotC)
  pplt.title('Density as function of delay')
  pplt.scatter(Delay*1E15, Ne_final, c=Efield1/1E9, label=r'$N_{exc}$[Delay]')
  # ax3.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
  pplt.xlabel('Double pulse delay (fs)')
  pplt.ylabel(r'Final density (m$^{-3}$)')
#  plt.ylim((0,5e28))
  ax3.set_yscale("log")

  #Plot 4: Irreversibility as function of delayD
  ax4 = plt.subplot(plotD)
  pplt.title('Irreversibility factor')
  pplt.scatter(Delay*1E15, irreversibility, c=Efield1/1E9, label='Irreversibility[Delay]')
  # ax4.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
  pplt.xlabel('Double pulse delay (fs)')
  pplt.ylabel('Irreversibility factor')
  pplt.ylim((0,1))

##TODO: shall we make a multi-wavelength figure?
##TODO: - If wavelength2!=wavelength1, yes.
##TODO: - If len(set(wavelength1)) > 1, yes.  

# Final saving
figmode=''

if(EnableDoublePulseAnalysis == 1):
  figmode = figmode+'-DoublePulse'
if(EnableInducedFieldAnalysis == 1):
  figmode = figmode+'-InducedFields'

figname = 'Metaplot'+figmode
plt.savefig(figname+'.eps', bbox_inches='tight')
plt.savefig(figname+'.png', bbox_inches='tight')
if(ShowFinalPlot==1):
  plt.show()


## Plot as function of wavelength for one given peak intensity
## Plot as function of peak intensity for one given wavelength.
## Plot as function of double pulse delay


