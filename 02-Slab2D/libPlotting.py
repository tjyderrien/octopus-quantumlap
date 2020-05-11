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

## @package libPlotting
# Define shared functions for plotting routines and database construction

from numpy import genfromtxt, loadtxt, chararray
import glob
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.cm as cm
from matplotlib import rc
from matplotlib.legend_handler import HandlerLine2D

from libString import *
from libAtomicUnits import *
import calculateVolume as Sample

Header = "[libPlotting]"

## Load geometry, n_el, laser and total_current files out of an Octopus TD run (atomic units)
# @param path: provides the link to the folder where we want to extract resuls from. 
def LoadRunOutputs(path):
  
  #geometry = path+"../Si.xyz" #static method
  #laser=loadtxt("../td.general/laser", dtype='str', skiprows=6) #should be replaced by genfromtxt method, as some errors are present in laser file
  
  #NOTE: the laser file must be treated by a script because there is a bug there.  
  try:
      laser=np.genfromtxt(path+"../td.general/laser.tmp", dtype='str', skip_header=6, loose=False, invalid_raise=True, filling_values=0.)
  except: 
      print Header+"** Critical error: Laser file could not be retrieved. Stopping the database fetching... Did you run plotAll.sh first? Laser file has some column width problem from octopus. "
      exit()
  
  # PostTreatedDensity = loadtxt("postCalculationNexcited.log", dtype='str', skiprows=1) #Unit: number of electrons, which should be multiplied by cell volume (AU). 
  TotalCurrent=loadtxt(path+"../td.general/total_current", dtype='str', skiprows=4)
  Energy = np.genfromtxt(path+"../td.general/energy", dtype='str', skip_header=5)
  try:
      # This should work for dielectrics
      Ne=loadtxt(path+"../td.general/n_ex", dtype='str', skiprows=9) #total number of excited electrons, not a density!
  except:
      # This should be for metallic materials.
      print Header+"** Capturing N_exc: failed --> maybe the material is a metal? Then we set N_exc=0 and set the rest using laser field."
      Ne=Energy #we just copy energy, to re-use the same routines for database construction. plotAll.sh will have to be adapted to avoid plotting n_ex. 
      Ne[:,2]=0. #we set the column with n_exc values equal to 0. We keep the timesteps of Ne, employed by other functions. 
  return Ne, laser, TotalCurrent, Energy

## Captures the value of a laser parameter in the Octopus input file located in FILEPATH.
# Hence, localize yourself in td.general before using this function. 
def GetParameterValueFromInputFile(filepath, string): #{{{
    Header="GetParameterValueFromInputFile: "
    wavelength_inp       = LookUpStringInFile(filepath, string)
    # print Header+str(wavelength_inp)
    # split it using = sign, take the rhs then split again. 
    wavelength_SI_tmp    = wavelength_inp.split("=")
    try:
        wavelength_SI_tmp2   = wavelength_SI_tmp[1] #rhs
    except: 
        print Header+"Selection of r.h.s. failed in INP. Error in:"+str(wavelength_SI_tmp)
        exit()
    wavelength_SI_tmp3   = wavelength_SI_tmp2.split(" ") #select the rhs, and split the result from spaces
    print wavelength_SI_tmp3
    try: 
      wavelength_SI   = float(wavelength_SI_tmp3[1])
    except: 
      wavelength_SI   = float(wavelength_SI_tmp3[2])
    return wavelength_SI
#}}}

## Lists the octopus runs available on HDD, and complete the existing database.  
# 0. Load the existing database. 
# 1. List all the available runs
# 2. For each run, check if it is already available in database by comparing the ID of the run. 
# 3. If not available, we shall add it to the database. Otherwise, we pass to the next one. 
# Keep care to not delete the already available database. 
#  @param path: specifies the paths to consider. Wilcards "*" are advised. Example: '*/*/td.general/'
def ListAllOutput(path='*/*/td.general/', DatabaseName='Database.csv', FromScratch=False): #{{{
  Header="[OctopusSlabs-libPlotting] "
  
  #Load the existing database                                                                    # 1. Load the existing database
  try: 
    Database = np.loadtxt(DatabaseName, skiprows=0)
    print Header+"** Info: the existing database was found. We will reuse the existing DB. "
    FromScratch = False
  except: 
    print Header+"** Error: the existing database is empty. We will generate new DB. "
    FromScratch = True
    Database = np.empty((0,32))

  PlottingTable = np.empty((0,32)) ##NOTE: how to define tau_pulse for 2 colors? 
  #PlottingTable = np.empty((0,15)) #missing TheoryLevel
  # List all output_iter files
  print Header+"Fetching all subfolders...\n"
  for FolderName in glob.glob(path): #{{{
    #TODO: this is not general enough, directory level does not help
    print Header+"** Info: entering folder: "+str(FolderName)+"\n"

    # Import the quantities from the subfolders
    try:
      Ne, laser, TotalCurrent, Energy = LoadRunOutputs(FolderName)
    except:
      print Header+"** Warning: One or more files could not be read. There are some problems on laser file in octopus. Execute plotAll.sh to fix it. Then come back to me. I jump this folder:"+str(FolderName)
      continue #could not load the files, it did not run, then jump this run.
    
    # Collect the laser parameters
    filename             = FolderName+"../inp"
    
    # We can start by checking the RUN ID. If not available, no need to register this run... 
    print Header+"** Capturing the RUN ID of octopus on Draco..."
    runid=LookUpLatestLogFile(FolderName+"../octopus-td")
    print Header+"LogFile="+str(runid)+", folder = "+FolderName 
    try: #there is quite much loss here, 30% maybe. 
      int(runid)
    except:
      print Header+"** Warning: JOB ID could not be captured. Jumping to next computation."
      continue
    
    TheoryLevel          = 1 #default value
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
    
    wavelength_SI        = GetParameterValueFromInputFile(filename, "lambdaA = ")
    PulseDuration_SI     = GetParameterValueFromInputFile(filename, "tauAsi = ")
    DoublePulseDelay_SI  = GetParameterValueFromInputFile(filename, "pulseDelaySI = ")
	 
    FieldInduced_inp     = LookUpStringInFile(filename, "GaugeFieldPropagate = ")
    GaugeVectorField_inp = LookUpStringReturnNextLine(filename, "GaugeVectorField")
    GaugeFieldDelay_inp  = LookUpStringInFile(filename, "GaugeFieldDelaySI =")
    
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
    try:
      CEP1_inp             = LookUpStringInFile(filename, "CEPA = ")
      print Header+"** CEP1_inp: "+str(CEP1_inp)
      CEP1_SI              = eval(CEP1_inp[7:]) #the wrong expression can be selected sometimes
      #if(len(CEP1_SI)==0): #there exists some cases where "try" succeeds but CEP1 is empty...
        #CEP1_SI = 0.5*pi
        #print Header+"** Exception on CEP1_SI."
    except: 
      CEP1_inp             = "pi/2."
      print Header+"** Exception: CEP1 was se to 0 due to failure."
      CEP1_SI              = 0.5*pi

    #CEP1_SI              = CEP1_SI - pi/2. #convention due to some shift in octopus #exec(CEP1_inp[7:])
    print Header+"** Info: CEP1: "+CEP1_inp[7:]
       
    #Efield2_inp          = LookUpStringInFile(filename, "EfieldBsi =")
    Efield1_SI_input     = GetParameterValueFromInputFile(filename, "EfieldAsi = ")
    Efield2_SI_input     = GetParameterValueFromInputFile(filename, "EfieldBsi = ") #If non-zero, we enable double-pulse analysis. Else: disabled. 
        
    if(Efield2_SI_input > 1.0):#{ #get the second pulse data
      print Header+"** Entering second pulse input."
      #wavelength2_inp       = LookUpStringInFile(filename, "lambdaB =")
      wavelength2_SI        = GetParameterValueFromInputFile(filename, "lambdaB = ")
	 
      #PulseDuration2_inp    = LookUpStringInFile(filename, "tauBsi = ")
      PulseDuration2_SI     = GetParameterValueFromInputFile(filename, "tauBsi = ")
      
      CEP2_inp              = LookUpStringInFile(filename, "CEPB =")
      CEP2_SI               = eval(CEP2_inp[7:])
      CEP2_SI               = CEP2_SI - pi/2. #convention due to some shift in octopus
      #}
    else: #{{{
      wavelength2_SI=wavelength_SI
      PulseDuration2_SI=PulseDuration_SI
      CEP2_SI=CEP1_SI
      print Header+"Info: Pulse 2 was not used. "
    #}}}
    
    ######## Identify the type of the pulse envelope 

    # Gaussian | SinSquared pulse | bugged Soft Top Hat pulse (2) | Fixed STH (3) ?
#    del PulseType_inp, PulseType_value
    PulseType_inp   = LookUpStringInFile(FolderName+"../inp", '| "FieldColorA-') #BUG: To which file do we refer each time ?
    PulseType_value = PulseType_inp[48:] #this allows to have unique solution
    print Header+"** Debug info: PulseType_value = "+str(PulseType_value)
    if("Gaussian" in PulseType_value):
      print Header+"** Info: Gaussian pulse detected..."
      PulseTypeName="gaussian"
    elif ("SinSquared" in PulseType_value):
      print Header+"** Info: Squared Sinus pulse detected..."
      PulseTypeName="SinSquared"
    elif ("SinFourth" in PulseType_value): 
      print Header+"** Info: Sinus**4 pulse detected..."
      PulseTypeName="SinFourth"
    elif ("STH" in PulseType_value):
      print Header+"** Info: BUGGY Soft Top Hat pulse detected..."
      PulseTypeName="STH"
    elif ("StH4" in PulseType_value): #fixed StH with sin4 envelope #comes first in priority
      print Header+"** Info: SIN4 soft top hat pulse detected..."
      PulseTypeName="StH4"
    elif ("StH" in PulseType_value): #fixed STH 
      print Header+"** Info: FIXED Soft Top Hat pulse detected..."
      PulseTypeName="StH"
    else:
      print Header+"** Error: pulse type was not recognized."
      print Header+"** Error PulseType_value = "+str(PulseType_value)
      exit()
    print Header+"** Info: PulseTypeID is fetching "+str(PulseTypeName)
    PulseType = PulseTypeID[str(PulseTypeName)] #refers to dictionnary of libPlotting
    del PulseType_inp, PulseType_value

    # Identify the RiseTime of the pulse (useful for STH and StH modes only)
    if(PulseType == 2): 
      RiseTime_position = FolderName.find("STH") #returns the position of STH chain in folder name
      RiseTime_name = FolderName[RiseTime_position+3:RiseTime_position+4] #have to capture it from the folder name right after "STH" but right before "/"
      try: 
        RiseTime = int(RiseTime_name)
        print Header+"** Info: Rise time of STH was "+str(RiseTime)+" optical cycle(s)."
      except: 
        RiseTime = 0
    elif(PulseType == 3): 
      RiseTime_position = FolderName.find("StH") #returns the position of STH chain in folder name
      RiseTime_name = FolderName[RiseTime_position+3:RiseTime_position+4] #have to capture it from the folder name right after "STH" but right before "/"
      try: 
        RiseTime = int(RiseTime_name)
        print Header+"** Info: Rise time of StH was "+str(RiseTime)+" optical cycle(s)."
      except: 
        RiseTime = 0
    else: 
      RiseTime = 0
    

    ######## Filter the induced fields #TODO: reshape these tests: lets filter in POST collection of data. Let just put 1 or 0 in the table. Then we can generate all the plots. 
    FieldInduced_enabled = 0 # default value
    if(FieldInduced_string == "#"):
      FieldInduced_enabled = 0
      print Header+"Info: InducedFields was disabled."
    elif(FieldInduced_string == "G" and FieldInduced_inp[22] == "y"):
      FieldInduced_enabled = 1
      print Header+"Info: InducedFields was enabled."
    else:
      FieldInduced_enabled = 0
      print Header+"Info: InducedFields was disabled."
    
    ### Checking if probe was enabled. Output the probe amplitude.
#    GaugeFieldAmplitude=0e0 #zero if it was not used.
    GaugeVectorField_tmp   = GaugeVectorField_inp.split('|') #structure: X|Y|Z #comment
    # if a "#" is found, then set 0!
    if "#" in GaugeVectorField_tmp[0]: 
        GaugeFieldAmplitude = 0e0
    else: 
        try: 
            GaugeFieldAmplitude = GaugeVectorField_tmp[0] #returns the X occurrence
        except:
            GaugeFieldAmplitude = 0e0
            print Header+"Warning: Selection of r.h.s. failed in INP. No data from "+str(GaugeVectorField_tmp)
    if(GaugeFieldAmplitude==''):
        del GaugeFieldAmplitude
        GaugeFieldAmplitude = 0e0
    
    GaugeFieldDelay=0e0 #zero if was not used
    GaugeFieldDelay_tmp    = GaugeFieldDelay_inp.split("=")
    try:
        GaugeFieldDelay    = GaugeFieldDelay_tmp[1] #rhs
    except: 
        GaugeFieldDelay=0e0
        print Header+"Warning: Selection of r.h.s. failed in INP. No data in:"+str(GaugeFieldDelay_tmp)
    
    ### Captures the functional which has been used for the computation.  #BUG: crashes when XCFunctional is not found. 
    Functional_value = "lda" #default value: lda_x + lda_c_pz_mod
    try: 
      Functional_string = LookUpStringInFile(filename, "XCFunctional")
    except:
      Functional_string = 'lda'
      print Header+"** Warning: search of XCFunctional failed. Default value was set: lda"
    
    if( len(Functional_string) > 0 ): #should be run only if Functional_string is defined. 
      Functional_value = Functional_string[13:] #capturing name of the functional
      Functional_lda_position = Functional_value.find("lda_x")
      Functional_mgga_position = Functional_value.find("MGGA")
      Functional_gga_position  = Functional_value.find("GGA")
      if(Functional_lda_position >= 0 and Functional_mgga_position == -1): #string was found
        Functional = 0 #lda
        print Header+"** Info: Functional LDA detected."
      elif(Functional_mgga_position >= 0 and Functional_lda_position == -1): #MMGA was enabled found
        Functional = 1 #mgga
        print Header+"** Info: Functional MGGA detected."
      elif(Functional_gga_position >= 0 and Functional_mgga_position == -1 and Functional_lda_position == -1): #GGA was found, not MGGA.
        Functional = 2 #GGA
        print Header+"** Info: Functional GGA detected."
      else: 
        print Header+"** Error: XCFunctional contained unknown case."+str(Functional_value)
        exit(-1)
    elif ( Functional_string == 'lda' ): #XCFunctional was not specified: assuming LDA.
      print Header+"Warning: XCFunctional was not specified. Assuming default value: lda."
      Functional = 0
    
    
    ### Captures if FreezeHXC was used. 
    FreezeHXC = 0
    try:
      FreezeHXC_string = LookUpStringInFile(filename, "TDFreezeHXC")
      print Header+"Info: TDFreezeHXC_string: "+str(FreezeHXC_string)
    except:
      print Header+"** Info: TDFreezeHXC was not found in the input file."
      #FreezeHXC_string = ''
    if("yes" in FreezeHXC_string): #We shall not make a strict test, which was failing. 
      FreezeHXC = 1
      print Header+"Info: TDFreezeHXC enabled."
    else:
      FreezeHXC = 0
      print Header+"Info: TDFreezeHXC not found."
    
    # Captures if MoveIon was enabled. 
    MoveIon = 0
    try: 
      MoveIonString = LookUpStringInFile(filename, "MoveIons")
      print Header+"Info: MoveIon_string: "+str(MoveIonString)
    except:
      print Header+"Info: MoveIon was not found in the input file."
    if("yes" in MoveIonString):
      MoveIon=1
    else:
      MoveIon=0
      print Header+"Info: MoveIon was not found in the input file."
    
    ### Remove contradictory case
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
    SampleCellVolume_AU = SuperCellVolume_AU #- Sample.VolumeOfVacuum(FolderName+"../") #Volume of the sample (without vacuum region) # WARNING: we disabled the SLAB volume computation !
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
    
    print Header+"** Importing simulation data..."
    
    # Catch the final excited electron density. 
    try: 
        timeNe           = np.asfarray(Ne[:,1]) #WARNING: this does exist for metals...
    except: 
        print Header+"** Error: It seems this run was not complete. We jump it."
        continue
    Nproj_AU         = np.asfarray(Ne[:,2]) #number of particles, not a density!
    density_proj_AU  = Nproj_AU / SampleCellVolume_AU #This is now a density (a.u.)
    density_proj_SI  = Density_AU_to_SI(density_proj_AU)
    sizeOfDensityTab = density_proj_SI.size

    density_proj_SI_final         = density_proj_SI[sizeOfDensityTab-1]
    density_proj_SI_max           = density_proj_SI.max()
    density_irreversibility_ratio = density_proj_SI_final / density_proj_SI_max
    
    print Header+"** Electron density analysis "
    print Header+"Final electron density = "+str(density_proj_SI_final)
    print Header+"Max.  electron density = "+str(density_proj_SI_max)
    print Header+"Irreversibility ratio (Final / Max) = "+str(density_irreversibility_ratio)
    ##==== Prepare plotting of energies ====
    timeEnergies   = np.asfarray(Energy[:,1]) #atomic units (! 1 is 2nd column)
    ElectronEnergy = np.asfarray(Energy[:,5]) #Hartree
    TotalEnergy    = np.asfarray(Energy[:,2]) #Hartree
    IonEnergy      = np.asfarray(Energy[:,3]) #Hartree
    Hartree        = np.asfarray(Energy[:,7]) #Hartree
    
    # Normalizing the database energies
    
    # We capture the energy increase due to the laser field. 
    # But there exist cases where the oscillations of the energy 
    # are larger than the energy amplitude. Therefore, lets just remove the ground state value. 
    ElectronEnergy = ElectronEnergy - ElectronEnergy[0]
    TotalEnergy    = TotalEnergy    - TotalEnergy[0]
    IonEnergy      = IonEnergy      - IonEnergy[0]
    Hartree        = Hartree        - Hartree[0]
    
    # Extracting the FINAL shift current (!!!)
    #timeCurrent          = np.asfarray(TotalCurrent[:, 1]) #atomic units
    TotalCurrent_X_AU     = np.asfarray(TotalCurrent[:,2]) #atomic units
    TotalCurrent_Y_AU     = np.asfarray(TotalCurrent[:,3]) #a. u.
    TotalCurrent_Z_AU     = np.asfarray(TotalCurrent[:,4]) #a. u. 
    TotalCurrent_size     = TotalCurrent_X_AU.size
    
    ## Selecting last value of the current
    TotalCurrent_X_AU_final = TotalCurrent_X_AU[TotalCurrent_size-1] 
    TotalCurrent_Y_AU_final = TotalCurrent_Y_AU[TotalCurrent_size-1]
    TotalCurrent_Z_AU_final = TotalCurrent_Z_AU[TotalCurrent_size-1]
    
    ## Also reporting maximum norm of the current
    TotalCurrent_AU_norm  = np.sqrt( np.power(TotalCurrent_X_AU,2)+np.power(TotalCurrent_Y_AU,2)+np.power(TotalCurrent_Z_AU,2) ) #computes norm of the current
    TotalCurrent_AU_max   = np.max(TotalCurrent_AU_norm)
    TotalCurrent_AU_Xmax  = np.max(TotalCurrent_X_AU)
    TotalCurrent_AU_Ymax  = np.max(TotalCurrent_Y_AU)
    TotalCurrent_AU_Zmax  = np.max(TotalCurrent_Z_AU)
    
    print Header+"******* WARNING: Here we neglect other polarization directions than X ! ******"
    print Header+"Final DC current: X: "+str(float(TotalCurrent_X_AU_final))+", Y: "+str(float(TotalCurrent_Y_AU_final))+", Z: "+str(float(TotalCurrent_Z_AU_final))+" a.u."
    TotalCurrent_SI_final = abs(TotalCurrent_X_AU_final) + abs(TotalCurrent_Y_AU_final) + abs(TotalCurrent_Z_AU_final)
    print Header+"** Missing: Final DC current: "+str(float(TotalCurrent_SI_final))+" TO BE CONVERTED TO S.I. "
    
    # Extracting the final energy increase (total, and electronic)
    # NOTE: Total energy is constant with time. Then we use total energy field for Hartree term for now. In future, we may add a field for kinetic energy of atoms. 
    HartreeEnergy_AU_final = Hartree[Hartree.size-1]
    HartreeEnergy_SI_final = float(Energy_Hartree_to_eV(HartreeEnergy_AU_final))
    
    ElectronEnergy_AU_final = ElectronEnergy[ElectronEnergy.size-1] 
    ElectronEnergy_SI_final = float(Energy_Hartree_to_eV(ElectronEnergy_AU_final))
    
    print Header+"Final electron energy (eV)     = "+str(ElectronEnergy_SI_final)
    print Header+"Final Hartree  energy (eV)     = "+str(ElectronEnergy_SI_final)
    # Append to the table
    # Structure: 
    # wavelength1, wavelength2, tau1, tau2, 
    # Epeak1, Epeak2, 
    # CEP1, CEP2, delay, 
    # kGrid, MeshSize, 
    # NeMax, NeEnd, density_irreversibility_ratio, EnergyMax
    # One line is one octopus run. O_o
     
    print Header+"Checking integrity before adding to table..."
    #print float(runid), float(wavelength_SI), float(PulseDuration_SI)
    #print float(Efield1_SI_input), float(CEP1_SI), float(wavelength2_SI)
    #print float(PulseDuration2_SI), float(Efield2_SI_input), float(CEP2_SI)
    #print float(DoublePulseDelay_SI), float(CrossPolarization_angle), float(TheoryLevel)
    #print float(FieldInduced_enabled)
    print "** Debug: GaugeFieldAmplitude: "+str(GaugeFieldAmplitude)
    #print float(GaugeFieldDelay)
    #print float(Functional),  float(FreezeHXC), float(MoveIon), float(density_proj_SI_max)
    #print float(density_proj_SI_final), float(density_irreversibility_ratio)
    #print float(TotalCurrent_SI_final), float(TotalCurrent_AU_max), float(TotalCurrent_AU_Xmax)
    #print float(TotalCurrent_AU_Ymax), float(TotalCurrent_AU_Zmax), float(FinalEnergy_SI_final)
    #print float(density_proj_SI_final/PulseDuration_SI), PulseType, RiseTime, float(timeNe[len(timeNe)-1])
	      
    StackToAdd = [ float(runid), float(wavelength_SI), float(PulseDuration_SI), 
	      float(Efield1_SI_input), float(CEP1_SI), float(wavelength2_SI), 
	      float(PulseDuration2_SI), float(Efield2_SI_input), float(CEP2_SI), 
	      float(DoublePulseDelay_SI), float(CrossPolarization_angle), float(TheoryLevel), 
	      float(FieldInduced_enabled), float(GaugeFieldAmplitude), float(GaugeFieldDelay), 
	      float(Functional),  float(FreezeHXC), float(MoveIon), 
	      float(density_proj_SI_max), float(density_proj_SI_final), float(density_irreversibility_ratio), 
	      float(TotalCurrent_SI_final), float(TotalCurrent_AU_max), float(TotalCurrent_AU_Xmax), 
	      float(TotalCurrent_AU_Ymax), float(TotalCurrent_AU_Zmax), float(ElectronEnergy_SI_final), 
	      float(HartreeEnergy_SI_final), float(density_proj_SI_final/PulseDuration_SI), PulseType, 
	      RiseTime, float(timeNe[len(timeNe)-1]) ]

    
    print Header+"Stack contains:"
    print StackToAdd
    print Header+"** Info: Sizeof(PlottingTable)="+str(np.shape(PlottingTable))+".\nSizeof(StackToAdd)="+str(np.shape(StackToAdd))
    
    if(FromScratch): #We construct a new database from scratch
      try:
        PlottingTable = np.vstack(( PlottingTable, StackToAdd )) 
        print Header+"** Info: these data were added to the stack."
      except:
        Header+"** Error: Failed to stack the database when looping on folders. Sizeof(PlottingTable)="+str(np.shape(PlottingTable))+". Sizeof(StackToAdd)="+str(np.shape(StackToAdd))
    else: #We only add the mising runs to the existing database
      # We have to check if this runs is not already available in already existing Database. 
      PreviousDB_ID = Database #Should be result from loadtxt. 
      #print PreviousDB_ID
      if int(runid) in PreviousDB_ID: #if this run was found in existing DB, we jump it. This formulation was successfully tested separately. 
          print Header+"** Info: ============================ THIS RUN WAS ALREADY EXISTING IN DB. Jumping it. ================="
          continue
      else: 
        PlottingTable = np.vstack(( PlottingTable, StackToAdd ))
        print Header+"** Info: these data were NEW and were ADDED to the stack."
    # PlottingTable should only contains FRESH data. 
    # Database contains the OLD data. 
    
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
  FreshData   = PlottingTable  
  NewDatabase = np.vstack(( Database, PlottingTable ))  
  return FreshData, NewDatabase
#}}} #end of function
  

## Unfold the data from an array of octopus results
# @param PlottingTable: array of floats containing the data to extract. 
# /!\ The size of the table must match the one given in ListAllOutput()
def ExtractData(PlottingTable): #{{{
  runid              =   PlottingTable[:,0]
  wavelengths1       =   PlottingTable[:,1]
  pulses1            =   PlottingTable[:,2]
  Efield1            =   PlottingTable[:,3]
  CEP1               =   PlottingTable[:,4]
  wavelengths2       =   PlottingTable[:,5]
  pulses2            =   PlottingTable[:,6]
  Efield2            =   PlottingTable[:,7]
  CEP2               =   PlottingTable[:,8]
  Delay              =   PlottingTable[:,9]
  CrossPolarization  =   PlottingTable[:,10]
  TheoryLevel        =   PlottingTable[:,11]
  FieldInduced       =   PlottingTable[:,12]
  GaugeFieldAmplitude=   PlottingTable[:,13]
  GaugeFieldDelay    =   PlottingTable[:,14]
  Functional         =   PlottingTable[:,15]
  FreezeHXC          =   PlottingTable[:,16]
  MoveIon            =   PlottingTable[:,17]
  Ne_max             =   PlottingTable[:,18]
  Ne_final           =   PlottingTable[:,19]
  irreversibility    =   PlottingTable[:,20]
  TotalCurrent_final =   PlottingTable[:,21]
  TotalCurrent_max   =   PlottingTable[:,22]
  TotalCurrent_AU_Xmax = PlottingTable[:,23]
  TotalCurrent_AU_Ymax = PlottingTable[:,24]
  TotalCurrent_AU_Zmax = PlottingTable[:,25]
  ElectronEnergy_max =   PlottingTable[:,26]
  TotalEnergy_max    =   PlottingTable[:,27]
  Wpi                =   PlottingTable[:,28]
  PulseType          =   PlottingTable[:,29]
  RiseTime           =   PlottingTable[:,30]
  FinalTime          =   PlottingTable[:,31]

  return runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, TotalEnergy_max, Wpi, PulseType, RiseTime, FinalTime
#}}}

IndexTable = {'runid': 0, 'wavelengths1': 1, 'pulses1': 2, 'Efield1': 3, 'CEP1': 4, 'wavelengths2': 5, 'pulses2': 6, 'Efield2': 7, 'CEP2': 8, 'Delay':  9, 'CrossPolarization': 10, 'TheoryLevel': 11, 'FieldInduced': 12, 'GaugeFieldAmplitude':13, 'GaugeFieldDelay':14, 'Functional': 15, 'FreezeHXC':16, 'MoveIon': 17, 'Ne_max': 18, 'Ne_final': 19, 'irreversibility': 20, 'Current_final': 21, 'Current_max': 22, 'CurrentX_max': 23, 'CurrentY_max': 24, 'CurrentZ_max': 25, 'ElectronEnergy_final': 26, 'Energy_max': 27, 'Wpi': 28, 'PulseType': 29, 'RiseTime': 30, 'FinalTime': 31}

## Counts the number of available laser parameters in a set of octopus results
# @param Database: array of octopus results. It can be filtered before using FilterDatabase functions. See libDatabase.py. 
def CountTheNumberOfLaserParameters(Database):#{{{
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, Energy_max, Wpi, PulseType, RiseTime, FinalTime = ExtractData(Database)
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
  Functional_uniq        = list(set(Functional))
  FreezeHXC_uniq         = list(set(FreezeHXC))
  MoveIon_uniq           = list(set(MoveIon))
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
  print "** Functionals             = "+str(len(Functional_uniq))
  print "** FreezeHXC               = "+str(len(FreezeHXC_uniq))
  print "** MoveIon                 = "+str(len(MoveIon_uniq))
  
  print "=== COUPLING EFFECTS ==="
  print "** Delays             size = "+str(len(Delay_uniq))
  print "** CrossPolarization  size = "+str(len(CrossPolarization_uniq))

  print "====== RESULTS ====="
  print "** Ne_final           size = "+str(len(Ne_final))
  print "** Irreversibility    size = "+str(len(irreversibility))
  print "** Total energy       size = "+str(len(Energy_max))
  
  return(wavelengths1_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses1_uniq, wavelengths2_uniq, CEP2_uniq, Efield2_uniq, pulses2_uniq)
#}}}


#Selects single-color pulses of various wavelengths, with same field amplitude from a database
# @param Database: the array with all the data, 
# @param Efield1 (V/m): desired field amplitude
# @param pulseDuration (s): desired pulse duration
# @param CEP: desired carrier envelope phase (pulse1)
# @param PulseType: 0: gaussian, 1: squared sin, 2: soft-top hat (bugged), 3: soft top hat (fixed)
# @param RiseTime: 0: not important. <n>: softened top hat with n cylces for increase and n for decrease.
def FilterSinglePulseDatabaseWithParameters_wavelengths(Database, Efield1, pulseDuration, CEP, PulseType=3, RiseTime=2, MoveIon=0): #{{{ 
        PlottingSinglePulses = Database  #get the full database
        PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
        PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, Efield1*0.99, IndexTable['Efield1'])
        PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, Efield1*1.01, IndexTable['Efield1'])
        PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
        PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0. 
        PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, CEP, IndexTable['CEP1']) #removes fields where CEP != 0. #NOTE: Can be done for Si-LDA. 
        PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['FreezeHXC']) #removes fields where FreezeHXC != 0. 
        PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['FieldInduced']) #removes fields where FreezeHXC != 0. 
        PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 1., IndexTable['TheoryLevel']) #removes runs where TheoryLevel != dft. 
        PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, MoveIon, IndexTable['MoveIon'])
        PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, PulseType, IndexTable['PulseType']) #keeps only pulses of type StH (3)
        PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, pulseDuration, IndexTable['pulses1']) #keeps only pulses of selected duration in argument
        if(RiseTime >= 1): #select by RiseTime only if it is defined
                PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, RiseTime, IndexTable['RiseTime'])
        # Inform the user of the contents.
        CountTheNumberOfLaserParameters(PlottingSinglePulses)
        return PlottingSinglePulses
#}}}

#Selects single-color pulses, with same field amplitude from a database
# @param Database: the array with all the data, 
# @param Efield1 (V/m): desired field amplitude
# @param pulseDuration (s): desired pulse duration
# @param CEP: desired carrier envelope phase (pulse1)
# @param PulseType: 0: gaussian, 1: squared sin, 2: soft-top hat (bugged), 3: soft top hat (fixed)
# @param RiseTime: 0: not important. <n>: softened top hat with n cylces for increase and n for decrease.
def FilterSinglePulseDatabaseWithParameters_CEP(Database, wavelength1, Efield1, pulseDuration, PulseType=3, RiseTime=2, MoveIon=0, InducedFields=0): #{{{ 
        PlottingSinglePulses = Database  #get the full database
        PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
        PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, Efield1*0.99, IndexTable['Efield1'])
        PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, Efield1*1.01, IndexTable['Efield1'])
        PlottingSinglePulses = FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
        PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, wavelength1, IndexTable['wavelengths1'])
        PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0. 
#        PlottingSinglePulses = FilterDatabaseGreaterThan(PlottingSinglePulses, 0., IndexTable['CEP1']) #removes fields where CEP < 0. Just for wrong buggy reasons. 
        PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 0., IndexTable['FreezeHXC']) #removes fields where FreezeHXC != 0. 
        PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, InducedFields, IndexTable['FieldInduced']) 
        PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, 1., IndexTable['TheoryLevel']) #removes runs where TheoryLevel != dft. 
        PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, MoveIon, IndexTable['MoveIon'])
        PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, PulseType, IndexTable['PulseType']) #keeps only pulses of type StH (3)
        PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, pulseDuration, IndexTable['pulses1']) #keeps only pulses of selected duration in argument
        if(RiseTime >= 1): #select by RiseTime only if it is defined
                PlottingSinglePulses = FilterDatabase(PlottingSinglePulses, RiseTime, IndexTable['RiseTime'])
        # Inform the user of the contents.
        CountTheNumberOfLaserParameters(PlottingSinglePulses)
        return PlottingSinglePulses
#}}}

# For historical reasons, 
# Gaussian is working but non-zero derivative in instant 0
# SinSquared is working and physical. 
# STH was not well controlled (stupid mistake)
# StH was well controlled. 
# StH4 is experimental.
PulseTypeID = {'gaussian': 0, 'SinSquared': 1, 'STH': 2, 'StH': 3, 'StH4': 4, 'SinFourth': 5}
PulseTypeID_rev = {0: 'gaussian', 1: 'SinSquared', 2: 'STH', 3: 'StH', 4: 'StH4', 5: 'SinFourth'}

## Splits the currents into several objects (all in atomic units)
def GetTotalCurrent(TotalCurrent): #{{{
#    print TotalCurrent[:,1]
    timeCurrent=np.asfarray(TotalCurrent[:,1]) #Unit a.u.
    # timeCurrent=timeCurrent[:timeCurrent.size-ShrinkBy]
#    print timeCurrent[2]
    CurrentX=np.asfarray(TotalCurrent[:,2])
#    print CurrentX[2]
    CurrentY=np.asfarray(TotalCurrent[:,3])
    CurrentZ=np.asfarray(TotalCurrent[:,4])
    return timeCurrent, CurrentX, CurrentY, CurrentZ
#}}}

## Splits the gauge field into several objects (all in atomic units)
def GetTotalGaugeField(obj): #{{{
    t, X, Y, Z = GetTotalCurrent(obj)
    return t, X, Y, Z
#}}}

## Damps a signal from FinalTimeDamping instant. 
# @param timeCurrent: table of instants (SI)
# @param CurrentX: signal to be processed
# @param wavelength_SI: fundamental wavelength of the signal
# @param FinalDampingInstant: instant (SI) where to apply the damping
def SignalDampR(timeCurrent, CurrentX, wavelength_SI, FinalDampingInstant): #{{{

    omega = 2.*np.pi*c/wavelength_SI
    DampingTime = 2.*np.pi/omega #tau_rise
    timeCurrent_SI = timeCurrent #* Time_AU_to_SI(1.)
    
    #FinalDampingInstant = tmax #timeCurrent_SI[maxSize/2] #this one might be more reasonable
    print Header+"** Info: HHG dampingTime: "+str(round(DampingTime*1E15))+" fs."
    print Header+"** Info: HHG final instant: "+str(round(FinalDampingInstant*1E15))+" fs."
    print Header+"** Info: omega = "+str(round(omega))

    # CurrentX_HHG = CurrentX * np.power(np.sin(0.5 * np.pi * ( timeCurrent - InitiateDampingInstant - DampingTime ) / DampingTime ), 2.) * heaviside(timeCurrent - InitiateDampingInstant, 0.5) * (1.-heaviside(timeCurrent - InitiateDampingInstant - DampingTime , 5.))
    damperB = np.sin(np.pi/2. * (timeCurrent_SI-FinalDampingInstant)/DampingTime)**2 * (np.heaviside(timeCurrent_SI-FinalDampingInstant+DampingTime,0.5)) * (1.-np.heaviside(timeCurrent_SI-FinalDampingInstant,0.5)) #useful for the [t2, FinalDampingInstant] interval. 
    damperA = 1.- np.heaviside(timeCurrent_SI-(FinalDampingInstant-DampingTime),0.5)
    CurrentX_HHG = np.multiply(CurrentX, (damperA + damperB))
    return CurrentX_HHG
#}}}

## Plots a heat-map color picture of randomly meshed data. Outputs the results on a regular grid for further usage. 
# @param X: 2D-matrix
# @param Y: 2D-matrix
# @param Z: f(X, Y)
# @param QuantityTitle: string with title (LaTeX accepted)
# @param filename:      string with filename (EPS and PNG will be generated)
def plot2dHeatMap(BrutX, BrutY, BrutZ, QuantityTitle, filename, ShowFinalPlot=True, PlotGradient=False): #{{{
  print len(BrutX)
  print len(BrutY)
  print len(BrutZ)
  X=list(set(BrutX))
  Y=list(BrutY)
  #Z=DensityMax2Dsorted
  #YY, XX = np.meshgrid(X, Y)
  # If the size is big enough, let's make a plot.
  #if(len(DensityMax2D) >= 2): #{{{
  if(len(BrutZ) >= 2): #{{{
      #try:
      triang = tri.Triangulation(BrutX, BrutY)
      triang.set_mask(np.hypot(BrutX[triang.triangles].mean(axis=1),
                               BrutY[triang.triangles].mean(axis=1))
                      < 0.15)
      # We can project the triangular mesh on a orthogonal by performing interpolation as done here. https://matplotlib.org/3.1.0/gallery/images_contours_and_fields/triinterp_demo.html
      refiner = tri.UniformTriRefiner(triang)
      tri_refi, z_test_refi = refiner.refine_field(BrutZ, subdiv=3)
      
      E1_t, E2_t = np.meshgrid( np.linspace(np.min(BrutX), np.max(BrutX), 100), np.linspace(np.min(BrutY), np.max(BrutY), 100) )
      #I1_t = .5*c*epsilon_0*np.power(E1_t*1E9,2.)
      #I2_t = .5*c*epsilon_0*np.power(E2_t*1E9,2.)
      interp_cubic = tri.CubicTriInterpolator(tri_refi, z_test_refi, kind='geom') #kind='min_E' | 'geom'
      Nexc_t     = interp_cubic(E1_t, E2_t) #Now we should have Nexc_t interpolated on a chosen E1, E2 grid. It could be energy also... Let's keep it generic. 
      Nexc_t_grad = interp_cubic.gradient(E1_t, E2_t) #should provide the slope of Nexc(E1, E2) in both directions E1 and E2 :-)
      
      ## Plotting the quantity...
      plt.figure()
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
      plt.triplot(triang, lw=0.1, color='grey')
      levels = np.linspace(np.min(BrutZ), np.max(BrutZ), 10)
      #cmap = cm.get_cmaps(name='terrain', lut=None)
      cmap = cm.get_cmap(name='Blues', lut=None)
      plt.tricontourf(tri_refi, z_test_refi, levels=levels, cmap=cmap)
      #plt.tricontour(tri_refi, z_test_refi, levels=levels,
              #colors=['0.25', '0.5', '0.5', '0.5', '0.5'],
              #linewidths=[1.0, 0.5, 0.5, 0.5, 0.5])
      plt.colorbar()
      plt.tight_layout()
        #except:
            #CS = -1
            #print Header+"** Warning: failed to produce one of the plots."
        
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
      
      if(PlotGradient): #{{{
        print Header+"** Info: Plotting the gradient..."
        print Header+"** Debug: Shape of I1: "+str(np.shape(E1_t))
        print Header+"** Debug: Shape of I2: "+str(np.shape(E2_t))
        print Header+"** Debug: Shape of the grad(interpolated field X): "+str(np.shape(Nexc_t_grad[0]))
        print Header+"** Debug: Shape of the grad(interpolated field Y): "+str(np.shape(Nexc_t_grad[1]))
        
        fig2=plt.figure()
        plt.title(r"$\nabla_{I_1}$"+QuantityTitle)
        plt.xlabel(r'$I_1$ (W/m$^{2}$)')
        plt.ylabel(r'$I_2$ (W/m$^{2}$)')
        CS = plt.contourf(E1_t, E2_t, Nexc_t_grad[0], cmap=plt.cm.Blues)
        plt.triplot(triang, lw=0.1, color='grey')
        plt.colorbar(CS)
        plt.tight_layout()
        plt.savefig(filename+'-GradientE1.eps')
        plt.savefig(filename+'-GradientE1.png')
        if(ShowFinalPlot==1):
            plt.show()
            
        fig3=plt.figure()
        plt.title(r"$\nabla_{I_2}$"+QuantityTitle)
        plt.xlabel(r'$I_1$ (W/m$^{2}$)')
        plt.ylabel(r'$I_2$ (W/m$^{2}$)')
        CS = plt.contourf(E1_t, E2_t, Nexc_t_grad[1], cmap=plt.cm.Blues)
        plt.triplot(triang, lw=0.1, color='grey')
        plt.colorbar(CS)
        plt.tight_layout()
        plt.savefig(filename+'-GradientE2.eps')
        plt.savefig(filename+'-GradientE2.png')
        if(ShowFinalPlot==1):
            plt.show()
        #}}}
      ## Extract the E1, E2, Nexc_t to prepares sections and post-treated results on a regular grid. Figure of (Energy_pulse1+2)/(Energy_pulse1+Energy_pulse2)
      return E1_t, E2_t, Nexc_t
  #}}}
#}}}


## Plots a heat-map color picture of randomly meshed data
# @param X: 2D-matrix
# @param Y: 2D-matrix
# @param Z: f(X, Y)
# @param QuantityTitle: string with title (LaTeX accepted)
# @param filename:      string with filename (EPS and PNG will be generated)
def plot2dHeatMap_HHG(BrutX, BrutY, BrutZ, QuantityTitle, filename, ShowFinalPlot=True, PlotGradient=False, PlotTriangles=True, xlabel='', ylabel='', cmap_ext='RdBu_r', number_of_levels=10, xlimit=[0,10], ylimit=[0,10], Egap=2.56*e): #{{{
#    cmap='Blues'
#    cmap='RdBu_r'
#    afmhot_r, hot, hot_r, bone_r, 
#    colormap="bone_r"

  print len(BrutX)
  print len(BrutY)
  print len(BrutZ)
  X=list(set(BrutX))
  Y=list(BrutY)

#  YY, XX = np.meshgrid(X, Y)
  # If the size is big enough, let's make a plot.
  #if(len(DensityMax2D) >= 2): #{{{
  if(len(BrutZ) >= 2): #{{{
      #try:
      triang = tri.Triangulation(BrutX, BrutY)
      triang.set_mask(np.hypot(BrutX[triang.triangles].mean(axis=1),
                               BrutY[triang.triangles].mean(axis=1))
                      < 0.15)
      # We can project the triangular mesh on a orthogonal by performing interpolation as done here. https://matplotlib.org/3.1.0/gallery/images_contours_and_fields/triinterp_demo.html
      refiner = tri.UniformTriRefiner(triang)
      tri_refi, z_test_refi = refiner.refine_field(BrutZ, subdiv=3)
      
      E1_t, E2_t = np.meshgrid( np.linspace(np.min(BrutX), np.max(BrutX), 100), np.linspace(np.min(BrutY), np.max(BrutY), 100) )
      #I1_t = .5*c*epsilon_0*np.power(E1_t*1E9,2.)
      #I2_t = .5*c*epsilon_0*np.power(E2_t*1E9,2.)
      interp_cubic = tri.CubicTriInterpolator(tri_refi, z_test_refi, kind='geom') #kind='min_E' | 'geom'
#      interp_linear= tri.LinearTriInterpolator(tri_refi, z_test_refi) #consumes much less. 
      Nexc_t     = interp_cubic(E1_t, E2_t) #Now we should have Nexc_t interpolated on a chosen E1, E2 grid. 
      Nexc_t_grad = interp_cubic.gradient(E1_t, E2_t) #should provide the slope of Nexc(E1, E2) in both directions E1 and E2 :-)
      
      ## Plotting the quantity...
      plt.figure()
      #plt.title(r'$N_{exc}(E_1,E_2)$, $\lambda_1=$'+str(wavelength1*1E9)+r'nm, $\lambda_2=$'+str(wavelength2*1E9)+'nm.')
      fig = plt.figure()
      #ax = fig.add_subplot(111, projection='3d')
      #plt.gca().set_aspect('equal')
      plt.title(QuantityTitle)
      if(xlabel==''): 
         plt.xlabel(r'$E_1$ (V/nm)')
      else: 
         plt.xlabel(xlabel)
      if(ylabel==""):
         plt.ylabel(r'$E_2$ (V/nm)')
      else: 
         plt.ylabel(ylabel)
      #CS = plt.imshow(XX, YY, DensityMax2Dsorted) #, cmap=plt.cm.Blues)
      #CS = plt.contourf(XX, YY, DensityMax2Dsorted, cmap=plt.cm.Blues)
      #CS = ax.scatter(BrutX, BrutY, BrutZ) #, cmap=plt.cm.Blues) #WORKS but limited to boring dots
      #CS = ax.plot_wireframe(XX, YY, DensityMax2Dsorted) #, cmap=plt.cm.Blues) #FAILS: WHY do we have ground values ?! 
      #CS = ax.plot_surface(XX, YY, DensityMax2DsortedShape) #, cmap=plt.cm.Blues) #FAILS: WHY do we have ground values ?! 
      #CS = ax.plot_trisurf(BrutX, BrutY, BrutZ, cmap=plt.cm.Blues, antialiased=True, edgecolor='none') #WORKS
      #ax.view_init(90, 90)
      #See https://matplotlib.org/gallery/images_contours_and_fields/tricontour_smooth_user.html#sphx-glr-gallery-images-contours-and-fields-tricontour-smooth-user-py
      if(PlotTriangles):
        plt.triplot(triang, lw=0.1, color='grey')
      levels = np.linspace(np.min(BrutZ), np.max(BrutZ), number_of_levels)
      #cmap = cm.get_cmaps(name='terrain', lut=None)
      cmap = cm.get_cmap(name=cmap_ext, lut=None)
      plt.tricontourf(tri_refi, z_test_refi, levels=levels, cmap=cmap)
      #plt.tricontour(tri_refi, z_test_refi, levels=levels,
              #colors=['0.25', '0.5', '0.5', '0.5', '0.5'],
              #linewidths=[1.0, 0.5, 0.5, 0.5, 0.5])
      plt.colorbar()
      plt.tight_layout()
        #except:
            #CS = -1
            #print Header+"** Warning: failed to produce one of the plots."
        
      #CS = ax.scatter(XXi, YYi, BicolorNeFinal) #, cmap=plt.cm.Blues) #WORKS
      
      #CS = ax.plot_trisurf(Xi, Yi, BicolorNeFinal, cmap=plt.cm.Spectral)
      
      #CS = ax.contour(XX, YY, Z, cmap=plt.cm.Blues) #WORKS
      #plt.colorbar(CS)
      
      # Limiting the view
      plt.xlim((xlimit[0], xlimit[1]))
      plt.ylim((ylimit[0], ylimit[1]))

      # Adding vertical bars for the gap energy
      for n in [1,2,3,4]:
        plt.axvline(x=Egap/e/n, linestyle=":")
        plt.axhline(y=1E6*h*c/(Egap/n), linestyle=':') 
        plt.axvline(x=n*Egap/e, linestyle=":")
        # plt.axhline(y=1E6*h*c/(Egap*n), linestyle=':') #not relevent

      # plt.savefig(filename+'.eps')
      plt.savefig(filename+'.png')
      print Header+"** Info: saved "+filename+".png."
      if(ShowFinalPlot==1):
        plt.show()
      
      if(PlotGradient):
        print Header+"** Info: Plotting the gradient..."
        print Header+"** Debug: Shape of I1: "+str(np.shape(E1_t))
        print Header+"** Debug: Shape of I2: "+str(np.shape(E2_t))
        print Header+"** Debug: Shape of the grad(interpolated field X): "+str(np.shape(Nexc_t_grad[0]))
        print Header+"** Debug: Shape of the grad(interpolated field Y): "+str(np.shape(Nexc_t_grad[1]))
        
        fig2=plt.figure()
        plt.title(r"$\nabla_{I_1}$"+QuantityTitle)
        plt.xlabel(r'$I_1$ (W/m$^{2}$)')
        plt.ylabel(r'$I_2$ (W/m$^{2}$)')
        CS = plt.contourf(E1_t, E2_t, Nexc_t_grad[0], cmap=plt.cm.Blues)
        plt.triplot(triang, lw=0.1, color='grey')
        plt.colorbar(CS)
        plt.tight_layout()
        plt.savefig(filename+'-GradientE1.eps')
        plt.savefig(filename+'-GradientE1.png')
        if(ShowFinalPlot==1):
            plt.show()
            
        fig3=plt.figure()
        plt.title(r"$\nabla_{I_2}$"+QuantityTitle)
        plt.xlabel(r'$I_1$ (W/m$^{2}$)')
        plt.ylabel(r'$I_2$ (W/m$^{2}$)')
        CS = plt.contourf(E1_t, E2_t, Nexc_t_grad[1], cmap=plt.cm.Blues)
        plt.triplot(triang, lw=0.1, color='grey')
        plt.colorbar(CS)
        plt.tight_layout()
        plt.savefig(filename+'-GradientE2.eps')
        plt.savefig(filename+'-GradientE2.png')
        if(ShowFinalPlot==1):
            plt.show()
#}}}
  
