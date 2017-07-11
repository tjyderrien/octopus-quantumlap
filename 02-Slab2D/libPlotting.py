#!/usr/bin/env python
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2017 T. J.-Y. Derrien
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


## @module libPlotting
## Define shared functions for plotting routines

from numpy import genfromtxt, loadtxt, chararray
import glob
import numpy as np

from libString import *
from libAtomicUnits import *
import calculateVolume as Sample

## Load geometry, n_el, laser and total_current files (atomic units)
def LoadRunOutputs(path):
  geometry = path+"../Si.xyz" #static method
  Ne=loadtxt(path+"../td.general/n_ex", dtype='str', skiprows=9) #total number of excited electrons, not a density!
  #laser=loadtxt("../td.general/laser", dtype='str', skiprows=6) #should be replaced by genfromtxt method, as some errors are present in laser file
  laser=np.genfromtxt(path+"../td.general/laser", dtype='str', skip_header=9, loose=False, invalid_raise=True, filling_values=0.)
  # PostTreatedDensity = loadtxt("postCalculationNexcited.log", dtype='str', skiprows=1) #Unit: number of electrons, which should be multiplied by cell volume (AU). 
  TotalCurrent=loadtxt(path+"../td.general/total_current", dtype='str', skiprows=4)
  Energy = np.genfromtxt(path+"../td.general/energy", dtype='str', skip_header=5)
  return geometry, Ne, laser, TotalCurrent, Energy

## Lists the available octopus runs, and build the table of results
#  @param path: specifies the paths to consider. Wilcards "*" are advised. Example: '*/*/output_iter/'
def ListAllOutput(path='*/*/output_iter/'): #{{{
  PlottingTable = np.empty((0,17))
  Header="[libPlotting] "
  #PlottingTable = np.empty((0,15)) #missing TheoryLevel
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
    
    ## Checking the theory level
    TheoryLevel_inp      = LookUpStringInFile(filename, "TheoryLevel =")
    try:
	 TheoryLevel_string   = TheoryLevel_inp[14:]
    except:
	 TheoryLevel_string   = ''
    if(TheoryLevel_string == "independent_particles"):
	 TheoryLevel = 0 #independent particles
    else: 
	 TheoryLevel = 1 #dft
	 
    ##TODO: Getting k-grid value
    
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
    print Header+"CrossPol"+str(CrossPolarization_angle)
    
    ## Import CEP data
    CEP1_inp             = LookUpStringInFile(filename, "CEPA =")
    print CEP1_inp
    CEP1_SI              = eval(CEP1_inp[7:]) 
    CEP1_SI              = CEP1_SI - pi/2. #convention due to some shift in octopus #exec(CEP1_inp[7:])
    
    Efield2_inp          = LookUpStringInFile(filename, "EfieldBsi =")
    Efield2_SI           = float(Efield2_inp[11:15]) #If non-zero, we enable double-pulse analysis. Else: disabled. 
        
    if(Efield2_SI > 1.0):#{ #get the second pulse data
	 print "** Entering second pulse input."
	 wavelength2_inp       = LookUpStringInFile(filename, "lambdaB =")
	 wavelength2_SI        = float(wavelength2_inp[10:17])
	 
	 PulseDuration2_inp    = LookUpStringInFile(filename, "tauBsi = ")
	 PulseDuration2_SI     = float(PulseDuration2_inp[8:15])
	 
	 CEP2_inp              = LookUpStringInFile(filename, "CEPB =")
	 CEP2_SI               = eval(CEP2_inp[7:])
	 CEP2_SI               = CEP2_SI - pi/2. #convention due to some shift in octopus
	 #}
    else: #{{{
	 wavelength2_SI=wavelength_SI
	 PulseDuration2_SI=PulseDuration_SI
	 CEP2_SI=CEP1_SI
    #}}}
    
    ######## Filter the induced fields #TODO: reshape these tests: lets filter in POST collection of data. Let just put 1 or 0 in the table. Then we can generate all the plots. 
    FieldInduced_enabled = 0
    if(FieldInduced_string == "#"):
	 FieldInduced_enabled = 0
    elif(FieldInduced_string == "G" and FieldInduced_inp[22] == "y"):
	 FieldInduced_enabled = 1
    else:
	 FieldInduced_enabled = 0
    ### Remove contradictory caes
    #if(EnableInducedFieldAnalysis == 1 and FieldInduced_string == "#"):
	 #print "Jumped: InducedFieldON but DISABLED in computation. [Case 1]"
	 #continue
    #if(EnableInducedFieldAnalysis == 0 and FieldInduced_string == "G" and FieldInduced_inp[22] == "y"):
	 #print "Jumped InducedFieldOFF but ENABLED in computation. [Case 2]"
	 #continue
    ### Matching cases also may have some exceptions
    #else: #{
	 #if(FieldInduced_string == "#"):
	   #print "Induced Field was DISABLED in computation."
	 #elif(FieldInduced_string == "G"):
	   #print "Induced field was ENABLED in computation."
	 #else: #{
	   #if(EnableInducedFieldAnalysis==1):
		#print "Induced field: contradiction. GaugeFieldPropagate was not found in this computation. Jumped."
		#continue
	   #else:
		#print "Induced field was not found in this computation. But we take it into account."
	##}
    ##}

    ## Double-pulse cases: jump selections which don't match with desire. 
    #if(Efield2_SI > 1E0 and EnableDoublePulseAnalysis == 0): #In case we found a double-pulse
	 #print "Jump [Case 3]"
	 #continue #we jump this run, then 
    #elif(Efield2_SI < 1E0 and EnableDoublePulseAnalysis == 1): 
	 #print "Jump [Case 4]"
	 #continue

  #  else: 
  #    print "SINGLE-PULSE MODE"
	 # Remove the pulses with double

    print Header+"Wavelength (nm): "+str(wavelength_SI)
    print Header+"Pulse duration (s): "+str(PulseDuration_SI)
    print Header+"Double pulse delay (s): "+str(DoublePulseDelay_SI)
    print Header+"Carrier Envelope Phase (deg): "+str(CEP1_SI * 180./pi)
    
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
    
    print Header+"Electric field amplitude 1 = "+str(Efield_amplitude_SI / 1E9)+" V/nm."
    print Header+"Electric field amplitude 2 = "+str(Efield2_amplitude_SI / 1E9)+" V/nm."
    
    # Capturing cell volume
    try:
	 SuperCellVolume_AU = Sample.VolumeOfSupercell_Log(FolderName+"../")
    except:
	 print Header+"** Warning: sample computation volume was jumped!"
	 continue
    SuperCellVolume_SI = Volume_AU_to_SI(SuperCellVolume_AU)
    SampleCellVolume_AU = SuperCellVolume_AU - Sample.VolumeOfVacuum(FolderName+"../") #Volume of the sample (without vacuum region) 
    CellVolume_AU = Sample.VolumePerPoint_Log(FolderName+"../") #0.06385 #0.07820 #Bohr**3
    SampleCellVolume_SI = Volume_AU_to_SI(SampleCellVolume_AU)
    
    print Header+"** Sample volumes **"
    print Header+"Super-cell total volume = "+str(SuperCellVolume_AU)+" Bohr**3."
    print Header+"Super-cell total volume = "+str(SuperCellVolume_SI)+" m^{3}." 
    print Header+""
    print Header+"Sample volume = "+str(SampleCellVolume_AU)+" Bohr**3."
    print Header+"Sample vacuum volume = "+str(SampleCellVolume_SI)+" m^{3}."
    print Header+""
    print Header+"Cell volume = "+str(CellVolume_AU)+" Bohr**3."
    print Header+"Cell volume = "+str(Volume_AU_to_SI(CellVolume_AU))+" m^{-3}."
    print Header+"" 
    
    print Header+"** Capturing the RUN ID of octopus on Draco..."
    runid=LookUpLatestLogFile(FolderName+"../")
    print Header+"LogFile="+str(runid)
    
    print Header+"** Importing simulation data..."
    
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
    
    TotalEnergy_size = ElectronEnergy.size #TotalEnergy.size #TODO: to clean
    TotalEnergy_AU_final = ElectronEnergy[TotalEnergy_size-1] # TotalEnergy[TotalEnergy_size-1] #selects last value of energy
    FinalEnergy_SI_final = float(Energy_Hartree_to_eV(TotalEnergy_AU_final))
    print "Total energy (eV)     = "+str(FinalEnergy_SI_final)
    # Append to the table
    # Structure: 
    # wavelength1, wavelength2, tau1, tau2, 
    # Epeak1, Epeak2, 
    # CEP1, CEP2, delay, 
    # kGrid, MeshSize, 
    # NeMax, NeEnd, density_irreversibility_ratio, EnergyMax
    # One line is one octopus run. O_o
     
    print "Checking integrity before adding to table..."
       
    StackToAdd = [float(runid), float(wavelength_SI), float(PulseDuration_SI), 
	      float(Efield_amplitude_SI), float(CEP1_SI), float(wavelength2_SI), 
	      float(PulseDuration2_SI), float(Efield2_amplitude_SI), float(CEP2_SI), 
	      float(DoublePulseDelay_SI), float(CrossPolarization_angle), float(TheoryLevel), 
	      float(FieldInduced_enabled), float(density_proj_SI_max), float(density_proj_SI_final), 
	      float(density_irreversibility_ratio), float(FinalEnergy_SI_final) ]
    
    print Header+"Stack contains:"
    print StackToAdd
    
    PlottingTable = np.vstack(( PlottingTable, StackToAdd ))
    
    # Cleaning the guys before next loop
    #del time, A1x, A2x, timestepTable, Efield, Efield2, geometry, Ne, laser, TotalCurrent, Energy
    #del wavelength1, wavelength2
    #del tau1, tau2
    #del Epeak1, Epeak2
    #del CEP1, CEP2
    #del delay, kGrid, MeshSize
    #del density_proj_SI_max, density_proj_SI_final, density_irreversibility_ratio
    #del EnergyMax
    #del timeEnergies, ElectronEnergy, TotalEnergy, IonEnergy, Hartree
    
    #}}}
  return PlottingTable
#}}} #end of function
  

## Unfold the data from an array of octopus results
# @param PlottingTable: array of floats containing the data to extract. 
# /!\ The size of the table must match the one given in ListAllOutput()
def ExtractData(PlottingTable): #{{{
  runid             =   PlottingTable[:,0]
  wavelengths1      =   PlottingTable[:,1]
  pulses1           =   PlottingTable[:,2]
  Efield1           =   PlottingTable[:,3]
  CEP1              =   PlottingTable[:,4]
  wavelengths2      =   PlottingTable[:,5]
  pulses2           =   PlottingTable[:,6]
  Efield2           =   PlottingTable[:,7]
  CEP2              =   PlottingTable[:,8]
  Delay             =   PlottingTable[:,9]
  CrossPolarization =   PlottingTable[:,10]
  TheoryLevel       =   PlottingTable[:,11]
  FieldInduced      =   PlottingTable[:,12]
  Ne_max            =   PlottingTable[:,13]
  Ne_final          =   PlottingTable[:,14]
  irreversibility   =   PlottingTable[:,15]
  Energy_max        =   PlottingTable[:,16]

  return runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Ne_max, Ne_final, irreversibility, Energy_max
#}}}

## Counts the number of available laser parameters in a set of octopus results
# @param Database: array of octopus results. It can be filtered before using FilterDatabase functions. See libDatabase.py. 
def CountTheNumberOfLaserParameters(Database):#{{{
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Ne_max, Ne_final, irreversibility, Energy_max = ExtractData(Database)
  # Generates lists of unique entries to count the range of parameters we have. 
  runid_uniq             = list(set(runid))
  wavelengths1_uniq      = list(set(wavelengths1))
  wavelengths2_uniq      = list(set(wavelengths2))
  Efield1_uniq           = list(set(Efield1))
  Efield2_uniq           = list(set(Efield2))
  pulses1_uniq           = list(set(pulses1))
  pulses2_uniq           = list(set(pulses2))
  Delay_uniq             = list(set(Delay))
  FieldInduced_uniq      = list(set(FieldInduced))
  TheoryLevel_uniq       = list(set(TheoryLevel))
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

  print "=== THEORY LEVEL ==="
  print "** Number of theory levels = "+str(len(TheoryLevel_uniq))
  print "** Induced Fields          = "+str(len(FieldInduced_uniq))
  
  print "=== COUPLING EFFECTS ==="
  print "** Delays             size = "+str(len(Delay_uniq))
  print "** CrossPolarization  size = "+str(len(CrossPolarization_uniq))

  print "====== RESULTS ====="
  print "** Ne_final           size = "+str(len(Ne_final))
  print "** Irreversibility    size = "+str(len(irreversibility))
  print "** Total energy       size = "+str(len(Energy_max))
  
  return(wavelengths1_uniq, CEP1_uniq, TheoryLevel_uniq)
#}}}

IndexTable = {'runid': 0, 'wavelengths1': 1, 'pulses1': 2, 'Efield1': 3, 'CEP1': 4, 'wavelengths2': 5, 'pulses2': 6, 'Efield2': 7, 'CEP2': 8, 'Delay':  9, 'CrossPolarization': 10, 'TheoryLevel': 11, 'FieldInduced': 12, 'Ne_max': 13, 'Ne_final': 14, 'irreversibility': 15, 'Energy_max': 16}

