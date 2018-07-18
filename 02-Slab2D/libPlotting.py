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


## @module libPlotting
## Define shared functions for plotting routines

from numpy import genfromtxt, loadtxt, chararray
import glob
import numpy as np

from libString import *
from libAtomicUnits import *
import calculateVolume as Sample

Header = "[libPlotting]"

## Load geometry, n_el, laser and total_current files out of an Octopus TD run (atomic units)
# @param path: provides the link to the folder where we want to extract resuls from. 
def LoadRunOutputs(path):
  geometry = path+"../Si.xyz" #static method
  Ne=loadtxt(path+"../td.general/n_ex", dtype='str', skiprows=9) #total number of excited electrons, not a density!
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
  return geometry, Ne, laser, TotalCurrent, Energy


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

  PlottingTable = np.empty((0,23)) ##NOTE: how to define tau_pulse for 2 colors? 
  #PlottingTable = np.empty((0,15)) #missing TheoryLevel
  # List all output_iter files
  print Header+"Fetching all subfolders...\n"
  for FolderName in glob.glob(path): #{{{
    #TODO: this is not general enough, directory level does not help
    print Header+"** Info: entering folder: "+str(FolderName)+"\n"

    # Import the quantities from the subfolders
    try:
      geometry, Ne, laser, TotalCurrent, Energy = LoadRunOutputs(FolderName)
    except:
      print Header+"** Warning: One or more files could not be read. There are some problems on laser file in octopus. Execute plotAll.sh to fix it. Then come back to me. I jump this folder:"+str(FolderName)
      continue #could not load the files, it did not run, then jump this run.
    
    # Collect the laser parameters
    filename             = FolderName+"../inp"
    
    # We can start by checking the RUN ID. If not available, no need to register this run... 
    print Header+"** Capturing the RUN ID of octopus on Draco..."
    runid=LookUpLatestLogFile(FolderName+"../")
    print Header+"LogFile="+str(runid)
    try:
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
       
    Efield2_inp          = LookUpStringInFile(filename, "EfieldBsi =")
    Efield2_SI           = float(Efield2_inp[11:15]) #If non-zero, we enable double-pulse analysis. Else: disabled. 
        
    if(Efield2_SI > 1.0):#{ #get the second pulse data
      print Header+"** Entering second pulse input."
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
    elif ("STH" in PulseType_value):
      print Header+"** Info: BUGGY Soft Top Hat pulse detected..."
      PulseTypeName="STH"
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
      Functional_mgga_position = Functional_value.find("mgga")
      if(Functional_lda_position >= 0 and Functional_mgga_position == -1): #string was found
        Functional = 0 #lda
        print Header+"** Info: Functional LDA detected."
      elif(Functional_mgga_position >= 0 and Functional_lda_position == -1): #MMGA was enabled found
        Functional = 1 #gga
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
    
    # Extracting the FINAL shift current (!!!)
    #timeCurrent          = np.asfarray(TotalCurrent[:, 1]) #atomic units
    TotalCurrent_X_AU     = np.asfarray(TotalCurrent[:,2]) #atomic units
    TotalCurrent_Y_AU     = np.asfarray(TotalCurrent[:,3]) #a. u.
    TotalCurrent_Z_AU     = np.asfarray(TotalCurrent[:,4]) #a. u. 
    TotalCurrent_size     = TotalCurrent_X_AU.size
    TotalCurrent_X_AU_final = TotalCurrent_X_AU[TotalCurrent_size-1] 
    TotalCurrent_Y_AU_final = TotalCurrent_Y_AU[TotalCurrent_size-1]
    TotalCurrent_Z_AU_final = TotalCurrent_Z_AU[TotalCurrent_size-1]
    print Header+"******* WARNING: Here we neglect other polarization directions than X ! ******"
    print Header+"Final DC current: X: "+str(float(TotalCurrent_X_AU_final))+", Y: "+str(float(TotalCurrent_Y_AU_final))+", Z: "+str(float(TotalCurrent_Z_AU_final))+" a.u."
    TotalCurrent_SI_final = abs(TotalCurrent_X_AU_final) + abs(TotalCurrent_Y_AU_final) + abs(TotalCurrent_Z_AU_final)
    print Header+"** Missing: Final DC current: "+str(float(TotalCurrent_SI_final))+" TO BE CONVERTED TO S.I. "
    
    # Extracting the total energy
    TotalEnergy_size = ElectronEnergy.size #TotalEnergy.size #TODO: to clean
    TotalEnergy_AU_final = ElectronEnergy[TotalEnergy_size-1] # TotalEnergy[TotalEnergy_size-1] #selects last value of energy
    FinalEnergy_SI_final = float(Energy_Hartree_to_eV(TotalEnergy_AU_final))
    print Header+"Total energy (eV)     = "+str(FinalEnergy_SI_final)
    # Append to the table
    # Structure: 
    # wavelength1, wavelength2, tau1, tau2, 
    # Epeak1, Epeak2, 
    # CEP1, CEP2, delay, 
    # kGrid, MeshSize, 
    # NeMax, NeEnd, density_irreversibility_ratio, EnergyMax
    # One line is one octopus run. O_o
     
    print Header+"Checking integrity before adding to table..."
       
    StackToAdd = [ float(runid), float(wavelength_SI), float(PulseDuration_SI), 
	      float(Efield_amplitude_SI), float(CEP1_SI), float(wavelength2_SI), 
	      float(PulseDuration2_SI), float(Efield2_amplitude_SI), float(CEP2_SI), 
	      float(DoublePulseDelay_SI), float(CrossPolarization_angle), float(TheoryLevel), 
	      float(FieldInduced_enabled), float(Functional),  float(FreezeHXC), 
	      float(density_proj_SI_max), float(density_proj_SI_final), float(density_irreversibility_ratio), 
	      float(TotalCurrent_SI_final), float(FinalEnergy_SI_final), float(0.5*density_proj_SI_final/PulseDuration_SI), PulseType, RiseTime ]
    
    print Header+"Stack contains:"
    print StackToAdd
    
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
  Functional         =   PlottingTable[:,13]
  FreezeHXC          =   PlottingTable[:,14]
  Ne_max             =   PlottingTable[:,15]
  Ne_final           =   PlottingTable[:,16]
  irreversibility    =   PlottingTable[:,17]
  TotalCurrent_final =   PlottingTable[:,18]
  Energy_max         =   PlottingTable[:,19]
  Wpi                =   PlottingTable[:,20]
  PulseType          =   PlottingTable[:,21]
  RiseTime           =   PlottingTable[:,22]

  return runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, Wpi, PulseType, RiseTime
#}}}

## Counts the number of available laser parameters in a set of octopus results
# @param Database: array of octopus results. It can be filtered before using FilterDatabase functions. See libDatabase.py. 
def CountTheNumberOfLaserParameters(Database):#{{{
  runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, Wpi, PulseType, RiseTime = ExtractData(Database)
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
  
  print "=== COUPLING EFFECTS ==="
  print "** Delays             size = "+str(len(Delay_uniq))
  print "** CrossPolarization  size = "+str(len(CrossPolarization_uniq))

  print "====== RESULTS ====="
  print "** Ne_final           size = "+str(len(Ne_final))
  print "** Irreversibility    size = "+str(len(irreversibility))
  print "** Total energy       size = "+str(len(Energy_max))
  
  return(wavelengths1_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses1_uniq, wavelengths2_uniq, CEP2_uniq, Efield2_uniq, pulses2_uniq)
#}}}

IndexTable = {'runid': 0, 'wavelengths1': 1, 'pulses1': 2, 'Efield1': 3, 'CEP1': 4, 'wavelengths2': 5, 'pulses2': 6, 'Efield2': 7, 'CEP2': 8, 'Delay':  9, 'CrossPolarization': 10, 'TheoryLevel': 11, 'FieldInduced': 12, 'Functional': 13, 'FreezeHXC':14, 'Ne_max': 15, 'Ne_final': 16, 'irreversibility': 17, 'Current_final': 18, 'Energy_max': 19, 'Wpi': 20, 'PulseType': 21, 'RiseTime':22}

# For historical reasons, 
# Gaussian is working but non-zero derivative in instant 0
# SinSquared is working and physical. 
# STH was not well controlled (stupid mistake)
# StH was well controlled. 
PulseTypeID = {'gaussian': 0, 'SinSquared': 1, 'STH': 2, 'StH': 3}

## Damps a signal from FinalTimeDamping instant. 
# @param timeCurrent: table of instants (at.u.)
# @param wavelength_SI: fundamental wavelength of the signal
# @param FinalDampingInstant: instant where to apply the damping
def SignalDampR(timeCurrent, wavelength_SI, FinalDampingInstant): #{{{

    omega = 2.*np.pi*c/wavelength_SI
    DampingTime = 2.*np.pi/omega #tau_rise
    timeCurrent_SI = timeCurrent * Time_AU_to_SI(1.)
    
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
