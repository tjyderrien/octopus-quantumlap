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

ShowFinalPlot              = 0

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
  


def ListAllOutput(path='*/*/output_iter/'): #{{{
  PlottingTable = np.empty((0,13))
  # List all output_iter files
  for FolderName in glob.glob(path): #{{{
    #TODO: this is not general enough, directory level does not help
    print FolderName

    # Import the quantities from the subfolders
    try:
	 geometry, Ne, laser, TotalCurrent, Energy = LoadRunOutputs(FolderName)
    except:
	 continue #could not load the files, so jump this run.
    
    # Collect the laser parameters
    filename             = FolderName+"../inp"
    
    wavelength_inp       = LookUpStringInFile(filename, "lambdaA =")
    wavelength_SI        = float(wavelength_inp[10:17])
    
    PulseDuration_inp    = LookUpStringInFile(filename, "tauAsi = ")
    PulseDuration_SI     = float(PulseDuration_inp[8:15])
	 
    DoublePulseDelay_inp = LookUpStringInFile(filename, "pulseDelaySI = ")
    DoublePulseDelay_SI  = float(DoublePulseDelay_inp[15:22])
    
    FieldInduced_inp     = LookUpStringInFile(filename, "GaugeFieldPropagate = ")
    try:
	 FieldInduced_string = FieldInduced_inp[0]
    except:
	 FieldInduced_string = ''
	 
    CrossPolarization_inp     = LookUpStringInFile(filename, "th2 =") #if not found, then empty
    if(CrossPolarization_inp):
	 CrossPolarization_angle = float(eval(CrossPolarization_inp[6:13]))
    else: 
	 CrossPolarization_angle = 0.
    print "CrossPol"+str(CrossPolarization_angle)
    
    ## Import CEP data
    CEP1_inp             = LookUpStringInFile(filename, "CEPA =")
    print CEP1_inp
    CEP1_SI              = float(pi/2.) #exec(CEP1_inp[7:])
    
    Efield2_inp          = LookUpStringInFile(filename, "EfieldBsi =")
    Efield2_SI           = float(Efield2_inp[11:15]) #If non-zero, we enable double-pulse analysis. Else: disabled. 
        
    if(Efield2_SI > 1.0):#{ #get the second pulse data
	 print "** Entering second pulse input."
	 wavelength2_inp       = LookUpStringInFile(filename, "lambdaB =")
	 wavelength2_SI        = float(wavelength2_inp[10:17])
	 
	 PulseDuration2_inp    = LookUpStringInFile(filename, "tauBsi = ")
	 PulseDuration2_SI     = float(PulseDuration2_inp[8:15])
	 
	 CEP2_inp              = LookUpStringInFile(filename, "CEPB =")
	 print CEP2_inp #       = exec(CEP2_inp[7:])
	 CEP2_SI               = float(pi/2.)
	 #}
    else: #{{{
	 wavelength2_SI=wavelength_SI
	 PulseDuration2_SI=PulseDuration_SI
	 CEP2_SI=CEP1_SI
    #}}}
    
    # Filter the induced fields #TODO: reshape these tests. 
    ## Remove contradictory caes
    if(EnableInducedFieldAnalysis == 1 and FieldInduced_string == "#"):
	 print "Jumped: InducedFieldON but DISABLED in computation. [Case 1]"
	 continue
    if(EnableInducedFieldAnalysis == 0 and FieldInduced_string == "G" and FieldInduced_inp[22] == "y"):
	 print "Jumped InducedFieldOFF but ENABLED in computation. [Case 2]"
	 continue
    ## Matching cases also may have some exceptions
    else: #{
	 if(FieldInduced_string == "#"):
	   print "Induced Field was DISABLED in computation."
	 elif(FieldInduced_string == "G"):
	   print "Induced field was ENABLED in computation."
	 else: #{
	   if(EnableInducedFieldAnalysis==1):
		print "Induced field: contradiction. GaugeFieldPropagate was not found in this computation. Jumped."
		continue
	   else:
		print "Induced field was not found in this computation. But we take it into account."
	#}
    #}

    # Double-pulse cases: jump selections which don't match with desire. 
    if(Efield2_SI > 1E0 and EnableDoublePulseAnalysis == 0): #In case we found a double-pulse
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
    print "Carrier Envelope Phase (deg): "+str(CEP1_SI * 180./pi)
    
    # Catch the peak field #We could also take the one of input file..
    time                = np.asfarray(laser[:,1]) #column 2
    A1x                 = np.asfarray(laser[:,2]) #column number 3
    A2x                 = np.asfarray(laser[:,5]) #column number 6
    timestepTable       = np.diff(time) #size N-1
    TimeStep            = timestepTable[0] #they are all equal, just select the first
    Efield              = - (np.diff(A1x)/np.diff(time)) #correct value is *1/c 
    Efield2             = - (np.diff(A2x)/np.diff(time))
    Efield_AU           = Efield/c_AU #fixing dimension #TODO: clean me that...
    Efield2_AU          = Efield2/c_AU
    Efield_amplitude_SI = np.max(abs(Field_AU_to_SI(Efield_AU)))
    Efield2_amplitude_SI= np.max(abs(Field_AU_to_SI(Efield2_AU)))
    
    print "Electric field amplitude 1 = "+str(Efield_amplitude_SI / 1E9)+" V/nm."
    print "Electric field amplitude 2 = "+str(Efield2_amplitude_SI / 1E9)+" V/nm."
    
    # Capturing cell volume
    try:
	 SuperCellVolume_AU = Sample.VolumeOfSupercell_Log(FolderName+"../")
    except:
	 continue
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
    timeNe           = np.asfarray(Ne[:,1])
    Nproj_AU         = np.asfarray(Ne[:,2]) #number of particles, not a density!
    density_proj_AU  = Nproj_AU / SampleCellVolume_AU #This is now a density (a.u.)
    density_proj_SI  = Density_AU_to_SI(density_proj_AU)
    sizeOfDensityTab = density_proj_SI.size

    density_proj_SI_final         = density_proj_SI[sizeOfDensityTab-1]
    density_proj_SI_max           = density_proj_SI.max()
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
    
    TotalEnergy_size = TotalEnergy.size
    TotalEnergy_AU_final = TotalEnergy[TotalEnergy_size-1] #selects last value of energy
    FinalEnergy_SI_final = Energy_Hartree_to_eV(TotalEnergy_AU_final)
    # Append to the table
    # Structure: 
    # wavelength1, wavelength2, tau1, tau2, 
    # Epeak1, Epeak2, 
    # CEP1, CEP2, delay, 
    # kGrid, MeshSize, 
    # NeMax, NeEnd, density_irreversibility_ratio, EnergyMax
    # One line is one octopus run. O_o

    # Simplified version first
    PlottingTable = np.vstack(( PlottingTable, [float(wavelength_SI), PulseDuration_SI, Efield_amplitude_SI, #0-2
				  CEP1_SI, float(wavelength2_SI), PulseDuration2_SI, #3-5
				  Efield2_amplitude_SI, CEP2_SI, DoublePulseDelay_SI, #6-8
				  CrossPolarization_angle, density_proj_SI_max, density_proj_SI_final, #11
				  density_irreversibility_ratio, 0.0 ] )) ##13
    
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
    
    #}}}
  return PlottingTable
  #}}} #end of function

# Capture the data
path='*/*/output_iter/'
PlottingTable = ListAllOutput(path)

# Save the table to a file
libDatabase.ExportToTxt(PlottingTable, 'Database.csv')

# Splitting the data
wavelengths1      = PlottingTable[:,0]
pulses1           = PlottingTable[:,1]
Efield1           = PlottingTable[:,2]
CEP1              = PlottingTable[:,3] #TODO: conflict of name
wavelengths2      = PlottingTable[:,4]
pulses2           = PlottingTable[:,5]
Efield2           = PlottingTable[:,6]
CEP2              = PlottingTable[:,7] #TODO: conflict of name
Delay             = PlottingTable[:,8]
CrossPolarization = PlottingTable[:,9]
Ne_max            = PlottingTable[:,10]
Ne_final          = PlottingTable[:,11]
irreversibility   = PlottingTable[:,12]
Energy_max        = PlottingTable[:,13]

#print wavelengths1.shape, pulses1.shape, Efield1.shape, Delay.shape, Ne_max.shape, Ne_final.shape, irreversibility.shape

print Header+"** Preparing the summary plots..."

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
print "** 3d plot: Efield1*           size = "+str(len(Efield1_uniq))
print "** 3d plot: wavelength1*       size = "+str(len(wavelengths1_uniq))
print "** 3d plot: pulse duration*    size = "+str(len(pulses1_uniq))
print "** 3d plot: CEP1               size = "+str(len(CEP1_uniq))

print "====== PULSE 2 ====="
print "** 3d plot: Efield2*           size = "+str(len(Efield2_uniq))
print "** 3d plot: wavelength2*       size = "+str(len(wavelengths2_uniq))
print "** 3d plot: pulse duration*    size = "+str(len(pulses2_uniq))
print "** 3d plot: CEP2               size = "+str(len(CEP2_uniq))

print "=== COUPLING EFFECTS ==="
print "** 3d plot: CrossPolarization  size = "+str(len(CrossPolarization_uniq))

print "====== RESULTS ====="
print "** 3d plot: Ne_final           size = "+str(len(Ne_final))
print "** 3d plot: Irreversibility    size = "+str(len(irreversibility))
print "** 3d plot: Total energy       size = "+str(len(Energy_max))

### Defining plotting structure

# Observables:        {N_exc^max, N_exc^final, irreversibility, Energy_max, Energy_final}
# Number of subplots: 

exit()

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
SizeY = 15

plt.figure(figsize=(SizeX,SizeY))

#Plot 1: Final density as function of delay
ax1 = plt.subplot(plotA)
if(EnableDoublePulseAnalysis==1): 
  # Number of sweeped parameters: 2
  plt.title('Final density (peak field, delay)')
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
  pplt.scatter(Efield1, Delay*1E15, s=34+30*Ne_final/Density, c=Ne_final/Density)
else:
  if(wavelength_numbers == 1): #{{{ #single color simulations
    pplt.title('Density as function of peak field')
    # Number of sweeped parameters: 3 (Efield, Delay, CrossPolarization)
    pplt.scatter(Efield1, Ne_final, c=Delay, s=(30.+CrossPolarization*50.), label=r'$N_{exc}$[Field]')
    # ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    pplt.xlabel('Field peak amplitude (V/m)')
    pplt.ylabel(r'Final density (m$^{-3}$)')
  #  pplt.ylim((1E16,4.*5e28))
    ax1.set_yscale('log')
  #}}}
  else: #{{{ #various runs of different colors
    pplt.title('Available simulation results')
    pplt.scatter(1E9*wavelengths1, Efield1/1E9, c=Efield1/1E9, label='List of runs')
    pplt.xlabel('Laser wavelength (nm)')
    pplt.ylabel('Electric field (V/m)')
#    pplt.ylim((0,1))
 
 #}}}

#Plot 2: irreversibility
ax2 = plt.subplot(plotB)

if(EnableDoublePulseAnalysis==1):
  plt.title('Irreversibility factor')
  pplt.ylabel('Delay (fs)')
#  pplt.y2label('Adibadicity parameter')
  pplt.xlabel('Field ampl. (V/m)')
#  area = np.po
  pplt.scatter(Efield1, Delay*1E15, s=34+32*irreversibility, c=irreversibility)
else:
  if(wavelength_numbers == 1): #{{{
    plt.title('Irreversibility factor')
    pplt.scatter(Efield1, irreversibility, c=Delay*1E15, label='Irreversibility')
    pplt.xlabel('Field peak amplitude (V/m)')
    pplt.ylabel('Irreversibility factor')
    pplt.ylim((0,1))
    ax22 = ax2.twinx()
    pplt.scatter(Efield1, gammaKeldysh(2.56*e, 0.228, Efield1, wavelengths1), c='r', label='Adiabaticity')
    pplt.ylabel('Adiabadicity factor')
    pplt.ylim((0,1))
    # ax2.set_yscale("log")
    ax2.legend(loc=2, borderaxespad=0.)
    # ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax22.legend(loc=1, borderaxespad=0.)
  #}}}
  else: #{{{
    pplt.title('Density as function of laser wavelength')
    pplt.scatter(1E9*wavelengths1, Ne_final, c=Efield1, s=(30.+CrossPolarization*50.), label=r'$N_{exc}$[Field]')
    pplt.xlabel('Laser wavelength (nm)')
    pplt.ylabel(r'Final density (m$^{-3}$)')
  #}}}
  
if(EnableDoublePulseAnalysis==1): #{{{
  if(wavelength_numbers == 1): 
    #Plot 3: Final density as function of delay
    pplt.title('Density as function of delay')
    pplt.scatter(Delay*1E15, Ne_final, s=5.+(Efield1/1E9)*10., c=10*np.rint(CrossPolarization*180./pi), label=r'$N_{exc}$[Delay], Pola='+str(CrossPolarization[0]*180./pi)+"deg")
    # ax3.legend(bbox_to_anchor=(0.00, -0.1), loc=3, borderaxespad=0.)
    pplt.xlabel('Double pulse delay (fs)')
    pplt.ylabel(r'Final density (m$^{-3}$)')
    # plt.ylim((1E16,4.*5e28))
    # plt.xlim((0.5E-1, 35.))
    # ax3.set_xscale("log")
    # ax3.set_yscale("log")

    #Plot 4: Irreversibility as function of delay
    ax4 = plt.subplot(plotD)
    pplt.title('Irreversibility factor')
    pplt.scatter(Delay*1E15, irreversibility, s=5.+(Efield1/1E9)*10., c=10*np.rint(CrossPolarization*180./pi), label='Irreversibility[Delay]')
    # ax4.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    pplt.xlabel('Double pulse delay (fs)')
    pplt.ylabel('Irreversibility factor')
    pplt.ylim((0,1))
  #}}}
else: #single pulse mode #TODO: this is redundant. Organization of various cases is sot clear. 
  if(wavelength_numbers > 1): #several colors
    pplt.scatter(1E9*wavelengths1, irreversibility, c=Efield1, label='Irreversibility')
    pplt.xlabel('Laser wavelength (nm)')
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


