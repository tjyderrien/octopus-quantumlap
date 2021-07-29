#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2016-2021 T. J.-Y. Derrien
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

import os.path
# import pdb
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from scipy.constants import pi, e, epsilon_0, c, h
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d, UnivariateSpline
import glob

import octopus_slabs.Samples.calculateVolume as  calcVol
import octopus_slabs.Libs.libAtomicUnits     as au
import octopus_slabs.Libs.libString          as libString
from   octopus_slabs.Libs.libLogging         import init_logger
from   spp_extended_theory.Libs              import libDatabase

Header = "[libPlotting]"

IndexFunctional={'LDA': 0, 'TB09': 3, 'GGA': 2, 'MGGA': 1}

logger = init_logger(__name__, verbose=True) #"plotFinalQuantities")

## Load geometry, n_el, laser and total_current files out of an Octopus TD run (atomic units)
# @param localpath: provides the link to the folder where we want to extract resuls from.
def LoadRunOutputs(localpath): #{{{
  
  #geometry = localpath+"../Si.xyz" #static method
  #laser=np.loadtxt("../td.general/laser", dtype='str', skiprows=6) #should be replaced by genfromtxt method, as some errors are present in laser file
  
  #NOTE: the laser file must be treated by a script because there is a bug there.  
  try:
      laser=np.genfromtxt(localpath+"../td.general/laser.tmp", dtype='str', skip_header=6, loose=False, invalid_raise=False, filling_values=0.)
  except:
      try: 
          laser = np.genfromtxt(localpath + "../td.general/laser", dtype='str', skip_header=6, loose=False,
                            invalid_raise=False, filling_values=0.)
          logger.warning(Header+"** Critical error: file could not be retrieved. Stopping the database fetching...\n"+
                             "Did you run plotAll.sh first? Laser file has some column width problem from octopus.\n"+
                            "We'll try to fly with usual laser file instead. ")
      except:
          logger.error(Header+"** Critical error: laser file could not be retrieved. It seems run is corrupted since\n"+
                              "laser nor laser.tmp could be loaded. \n"+
                              "This RUN should be DISABLED, maybe.") 
  
  # PostTreatedDensity = np.loadtxt("postCalculationNexcited.log", dtype='str', skiprows=1) #Unit: number of electrons, which should be multiplied by cell volume (AU).
  TotalCurrent=np.loadtxt(localpath+"../td.general/total_current", dtype='str', skiprows=4)
  Energy = np.genfromtxt(localpath+"../td.general/energy", dtype='str', skip_header=5)
  try:
      # This should work for dielectrics
      Ne=np.loadtxt(localpath+"../td.general/n_ex", dtype='str', skiprows=9) #total number of excited electrons, not a density!
  except:
      # This should be for metallic materials.
      logger.error((Header+"** Capturing N_exc: failed --> maybe the material is a metal? Then we set N_exc=0 and set the rest using laser field."))
      Ne=Energy #we just copy energy, to re-use the same routines for database construction. plotAll.sh will have to be adapted to avoid plotting n_ex. 
      Ne[:,2]=0. #we set the column with n_exc values equal to 0. We keep the timesteps of Ne, employed by other functions. 
  return Ne, laser, TotalCurrent, Energy
#}}}

## Captures the value of a laser parameter in the Octopus input file located in FILEPATH.
# Hence, localize yourself in td.general before using this function. 
def GetParameterValueFromInputFile(filepath, string): #{{{
    Header="GetParameterValueFromInputFile: "
    wavelength_inp       = libString.LookUpStringInFile(filepath, string)
    # print Header+str(wavelength_inp)
    # split it using = sign, take the rhs then split again. 
    wavelength_SI_tmp    = wavelength_inp.split("=")
    try:
        wavelength_SI_tmp2   = wavelength_SI_tmp[1] #rhs
    except: 
        logger.error((Header+"Selection of r.h.s. failed in INP. Error in:"+str(wavelength_SI_tmp)))
        exit()
    wavelength_SI_tmp3   = wavelength_SI_tmp2.split(" ") #select the rhs, and split the result from spaces
    logger.debug(wavelength_SI_tmp3)
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
#  @param localpath: specifies the paths to consider. Wilcards "*" are advised. Example: '*/*/td.general/'
def ListAllOutput(localpath='*/*/*/td.general/', DatabaseName='Database.csv', FromScratch=False): #{{{
  Header="[OctopusSlabs-libPlotting] "

  def GetMachineName(file):
      machinename_t = libString.LookUpStringInFile(file, "swim")
      logger.debug((Header + "** Debug: " + str(machinename_t)))
      try:
          MachineName = machinename_t.split()[5]
          logger.info((Header + " ** Name of the machine which produced the result: " + str(MachineName)))
      except:
          logger.warning((Header + "** Warning: machine name could not be extracted."))
      # MachineNameList={0: "EOS", 1: "Draco", 2: "Salomon", 3: "Prometheus"}
      # MachineNameHeader={"draco": 1, "salomon": 2, "p": 3}
      # if(MachineName)
      try:
          RunDate_t = libString.LookUpStringInFile(file, "started")
          RunDate = float(RunDate_t.strip().split()[3].replace("/", ""))
          logger.info((Header + "** Date of the run: " + str(RunDate)))
      except:
          RunDate = float("nan")
      return RunDate

  def GetTheoryLevel(filename):
      TheoryLevel = 1  # default value
      ## Checking the theory level
      TheoryLevel_inp = libString.LookUpStringInFile(filename, "TheoryLevel =")
      try:
          TheoryLevel_string = TheoryLevel_inp[14:]
      except:
          TheoryLevel_string = ''
      if (TheoryLevel_string == "independent_particles"):
          TheoryLevel = 0  # independent particles
      else:
          TheoryLevel = 1  # dft
      return TheoryLevel

      # Capturing real-mesh grid spacing
  def GetRealMeshGrid(filename):
      realmesh_inp = libString.LookUpStringReturnNextLine(filename, "%Spacing")
      realmesh_inp2 = realmesh_inp.strip()
      realmesh_inp3 = realmesh_inp2.split("|")
      realmesh_inp4 = realmesh_inp3[0]
      realmeshgrid = float(realmesh_inp4)
      # del realmesh_inp, realmesh_inp2, realmesh_inp3, realmesh_inp4
      return realmeshgrid

  ## Capture angle of secondPulse polarization in respect of first pulse
  def GetCrossPolarizationAngle(filename):
      # === Capture of the pulse2-angle of field
      CrossPolarization_inp = libString.LookUpStringInFile(filename, "th2 =")  # if not found, then zero
      # Since Dec 2020, one consider angle th1 for orienting the electric field in the crystal. In some inp files, th2 is set equal to th1. The following lines should account for this case.
      if (CrossPolarization_inp):
          if (type(CrossPolarization_inp[6:13]) == "float"):
              CrossPolarization_angle = float(eval(CrossPolarization_inp[6:13].strip()))
          elif ("th1" in CrossPolarization_inp[6:13]):  # if it is refering to th1
              CrossPolarization_angle = polarization_angle
          else:
              CrossPolarization_angle = 0.
      else:  # if th2 was not found in the input file
          CrossPolarization_angle = 0.
      logger.debug((Header + "CrossPol" + str(CrossPolarization_angle)))
      return CrossPolarization_angle

  ## Import CEP
  def GetCEP1(filename):
      try:
          CEP1_inp = libString.LookUpStringInFile(filename, "CEPA = ")
          logger.info((Header + "** CEP1_inp: " + str(CEP1_inp)))
          CEP1_SI = eval(CEP1_inp[7:])  # the wrong expression can be selected sometimes
          # if(len(CEP1_SI)==0): #there exists some cases where "try" succeeds but CEP1 is empty...
          # CEP1_SI = 0.5*pi
          # print Header+"** Exception on CEP1_SI."
      except:
          CEP1_inp = "pi/2."
          logger.warning((Header + "** Exception: CEP1 was se to 0 due to failure."))
          CEP1_SI = 0.5 * pi

      # CEP1_SI              = CEP1_SI - pi/2. #convention due to some shift in octopus #exec(CEP1_inp[7:])
      logger.info((Header + "** Info: CEP1: " + CEP1_inp[7:]))

      return CEP1_SI

  #Load the existing database
  # 1. Load the existing database
  try: 
    Database = np.loadtxt(DatabaseName, skiprows=0)
    logger.info((Header+"** Info: the existing database was found. We will reuse the existing DB. "))
    FromScratch = False
  except: 
    logger.warning((Header+"** Warning: the existing database is empty. We will generate new DB. "))
    FromScratch = True
    Database = np.empty((0,47))

  PlottingTable = np.empty((0,47)) ##NOTE: how to define tau_pulse for 2 colors?
  #PlottingTable = np.empty((0,15)) #missing TheoryLevel
  # List all output_iter files
  logger.info((Header+"Fetching all subfolders...\n"))
  for FolderName in glob.glob(localpath): #{{{
    #TODO: this is not general enough, directory level does not help
    logger.info((Header+"** Info: entering folder: "+str(FolderName)+"\n"))
    
    ## HERE, if we find a file named "DISABLED" in the folder, then we shall jump the folder, to avoid polluting the database. 
    # print(glob.glob(path+"../DISABLED")) #! Finds ALL occurrences of DISABLED files 
    if(os.path.isfile(FolderName+"../DISABLED")):
        logger.warning((Header+" ** WARNING: THIS RUN WAS REMOVED after finding DISABLED file in the folder. This prevents this run to pollute the database with, e.g., unconverged data. Note that this DISABLED file may have been added manually for good reasons. Good practice may be to write inside this file the reason for having disabled this run. "))
        continue

    # Import the quantities from the subfolders
    try:
      Ne, laser, TotalCurrent, Energy = LoadRunOutputs(FolderName)
    except:
      logger.warning((Header+"** Warning: One or more files could not be read. There are some problems on laser file in octopus. Execute plotAll.sh to fix it. Then come back to me. I jump this folder: "+str(FolderName)))
      continue #could not load the files, it did not run, then jump this run.
    
    # Collect the laser parameters in the input file name
    filename             = FolderName+"../inp"
    
    # We can start by checking the RUN ID. If not available, no need to register this run... 
    logger.info((Header+"** Capturing the RUN ID of octopus..."))
    runid=libString.LookUpLatestLogFile(FolderName+"../octopus-td")
    logger.info((Header+"LogFile="+str(runid)+", folder = "+FolderName))
    try: #there is quite much loss here, 30% maybe. 
      int(runid)
    except:
      logger.warning((Header+"** Warning: JOB ID could not be captured. Jumping to next computation."))
      continue
    
    LOG_newest = max(glob.iglob(FolderName+'../octopus*.log'), key=os.path.getmtime)

    # Capturing metadata
    RunDate = GetMachineName(LOG_newest)
    TheoryLevel = GetTheoryLevel(filename)
    kpointgrid=libString.CaptureKPointGridSizeFromINP(filename)
    realmeshgrid = GetRealMeshGrid(filename)

    # Capturing basic parameters
    wavelength_SI        = GetParameterValueFromInputFile(filename, "lambdaA = ")
    PulseDuration_SI     = GetParameterValueFromInputFile(filename, "tauAsi = ")
    DoublePulseDelay_SI  = GetParameterValueFromInputFile(filename, "pulseDelaySI = ")

    FieldInduced_inp     = libString.LookUpStringInFile(filename, "GaugeFieldPropagate = ")
    GaugeVectorField_inp2= libString.LookUpStringInFile(filename, "GaugeVectorField")
    GaugeVectorField_inp = libString.LookUpStringReturnNextLine(filename, "GaugeVectorField")
    GaugeFieldDelay_inp  = libString.LookUpStringInFile(filename, "GaugeFieldDelaySI =")
    
#    try:
#        FieldInduced_string = FieldInduced_inp[0]
#    except:
#        FieldInduced_string = ''

    #=== Capture the angle of the field
    polarization_angle_inp = libString.LookUpStringInFile(filename, "th1 =")
    if(polarization_angle_inp): #if found
         polarization_angle = float(eval(polarization_angle_inp[6:13]))
    else:
         polarization_angle = 0. #float('nan'). 0 for now.

    CrossPolarization_angle = GetCrossPolarizationAngle(filename) #angle between the 2 fields!
    
    CEP1_SI=GetCEP1(filename) #WARNING: GetCEP1 does not remove pi/2 but CEP2 does?
       
    #Efield2_inp          = libString.LookUpStringInFile(filename, "EfieldBsi =")
    Efield1_SI_input     = GetParameterValueFromInputFile(filename, "EfieldAsi = ")
    Efield2_SI_input     = GetParameterValueFromInputFile(filename, "EfieldBsi = ") #If non-zero, we enable double-pulse analysis. Else: disabled. 
        
    if(Efield2_SI_input > 1.0):#{ #get the second pulse data
      logger.info((Header+"** Entering second pulse input."))
      #wavelength2_inp       = libString.LookUpStringInFile(filename, "lambdaB =")
      wavelength2_SI        = GetParameterValueFromInputFile(filename, "lambdaB = ")
         
      #PulseDuration2_inp    = libString.LookUpStringInFile(filename, "tauBsi = ")
      PulseDuration2_SI     = GetParameterValueFromInputFile(filename, "tauBsi = ")
      
      CEP2_inp              = libString.LookUpStringInFile(filename, "CEPB =")
      CEP2_SI               = eval(CEP2_inp[7:])
      CEP2_SI               = CEP2_SI - pi/2. #convention due to some shift in octopus
      #}
    else: #{{{
      wavelength2_SI=wavelength_SI
      PulseDuration2_SI=PulseDuration_SI #NOTE: would be useful to define 0 here, if no second pulse was used. 
      CEP2_SI=CEP1_SI
      logger.info((Header+"Info: Pulse 2 was not used. "))
    #}}}
    
    ######## Identify the type of the pulse envelope 

    # Gaussian | SinSquared pulse | bugged Soft Top Hat pulse (2) | Fixed STH (3) ?
#    del PulseType_inp, PulseType_value
    PulseType_inp   = libString.LookUpStringInFile(FolderName+"../inp", '| "FieldColorA-')
    PulseType_value = PulseType_inp[48:] #this allows to have unique solution although its really rough
    logger.debug((Header+"** Debug info: PulseType_value = "+str(PulseType_value)))


    if("Gaussian" in PulseType_value):
        logger.info((Header+"** Info: Gaussian pulse detected..."))
        PulseTypeName="gaussian"
    elif ("SinSquared" in PulseType_value):
        logger.info((Header+"** Info: Squared Sinus pulse detected..."))
        PulseTypeName="SinSquared"
    elif ("SinFourth" in PulseType_value):
        logger.info((Header+"** Info: Sinus**4 pulse detected..."))
        PulseTypeName="SinFourth"
    elif ("STH" in PulseType_value):
        logger.info((Header+"** Info: BUGGY Soft Top Hat pulse detected..."))
        PulseTypeName="STH"
    elif ("StH4" in PulseType_value): #fixed StH with sin4 envelope #comes first in priority
        logger.info((Header+"** Info: SIN4 soft top hat pulse detected..."))
        PulseTypeName="StH4"
    elif ("StH" in PulseType_value): #fixed STH
        logger.info((Header+"** Info: FIXED Soft Top Hat pulse detected..."))
        PulseTypeName="StH"
    else:
        logger.error((Header+"** Error: pulse type was not recognized."))
        logger.error((Header+"** Error PulseType_value = "+str(PulseType_value)))
        exit()
    # TODO: test the match developed below, that uses PulseTypeID_rev
    # for element in PulseTypeID_rev: #for eah type of pulse we know
    #     if(element in PulseType_value):
    #         PulseTypeName = element
    #         logger.info(Header+"** Info: the pulse type is: "+str(element))
    logger.info(Header+"** Info: PulseTypeID is fetching "+str(PulseTypeName))
    PulseType = PulseTypeID[str(PulseTypeName)] #refers to dictionnary of libPlotting
    del PulseType_inp, PulseType_value

    # Identify the RiseTime of the pulse (useful for STH and StH modes only)
    if(PulseType == 2): 
      RiseTime_position = FolderName.find("STH") #returns the position of STH chain in folder name
      RiseTime_name = FolderName[RiseTime_position+3:RiseTime_position+4] #have to capture it from the folder name right after "STH" but right before "/"
      try: 
        RiseTime = int(RiseTime_name) #BUG: this method is no longer working. Where to get the RiseTime? Maybe in tauRise = 2*pi/omega in inp ?
        logger.info((Header+"** Info: Rise time of STH was "+str(RiseTime)+" optical cycle(s)."))
      except: 
        RiseTime = 0
    elif(PulseType == 3): 
      RiseTime_position = FolderName.find("StH") #returns the position of STH chain in folder name
      RiseTime_name = FolderName[RiseTime_position+3:RiseTime_position+4] #have to capture it from the folder name right after "STH" but right before "/" #BUG: this method is not working for the recent runs. It works with the old times when we had tested the RiseTime. 
      try: 
        RiseTime = int(RiseTime_name)
        logger.info((Header+"** Info: Rise time of StH was "+str(RiseTime)+" optical cycle(s)."))
      except: 
        RiseTime = 2
    else: 
      RiseTime = 0
    

    ######## Filter the induced fields #TODO: reshape these tests: lets filter in POST collection of data. Let just put 1 or 0 in the table. Then we can generate all the plots. 

## BOOLEAN: GaugeFieldPropagate <-- GaugeFieldPropagate (FieldInduced_enabled) OR GaugeVectorField

## 1. Define the simple boolean from input file (intent) given by GaugeFieldPropagate of the Input file.
    FieldInduced_enabled = 0 # default value
#    if("#" in FieldInduced_string[0:5]):
    if(len(FieldInduced_inp)>0):
        FieldInduced_inp_tab = FieldInduced_inp.split(" ")
        if("#" in FieldInduced_inp[0:5]): #Induced fields commented?
            FieldInduced_enabled = 0
            logger.info((Header+"Info: InducedFields was disabled."))
        elif("yes" in FieldInduced_inp_tab[2]): # "Yes is found right after sign "="
            FieldInduced_enabled = 1
            logger.info((Header+"Info: InducedFields was enabled [FieldInduced_inp="+str(FieldInduced_inp_tab)+"]."))
    #    else: #any other case should check if GaugeVectorField is present
    #      FieldInduced_enabled = 0
    #      print Header+"Info: InducedFields was disabled."
      
## 2. Then define a simple boolean checking if GaugeVectorField is present. WARNING: even if amplitude is zero, octopus will enable GaugeFieldPropagate.

    ### Checking if probe was enabled. Output the probe amplitude.
    GaugeVectorField_enabled=0
    if "#" in GaugeVectorField_inp2[0:5]: #NOTE: GaugeVectorField_inp2 contains the keyword declaration  (not the amplitude). Then GaugeFieldAmplitude is disabled
        GaugeVectorField_enabled = 0
        logger.info((Header+"** Info: GaugeVectorField was commented [GaugeVectorField_inp2="+str(GaugeVectorField_inp2)+"]."))
        GaugeFieldAmplitude=0e0 #zero if GaugeVectorField is undefined (useful for database)
    elif("GaugeVectorField" in GaugeVectorField_inp2): #there is no # in line beginning and GaugeVectorField is ON
        GaugeVectorField_enabled = 1
        logger.info((Header+"** Info: GaugeVectorField was enabled [GaugeVectorField_inp2="+str(GaugeVectorField_inp2)+"."))
    
## 3. Now we re-affect FieldInduced_enabled as follows. 
    FieldInduced_enabled = (GaugeVectorField_enabled or FieldInduced_enabled)
    del GaugeVectorField_enabled
    GaugeVectorField_tmp   = GaugeVectorField_inp.split('|') #structure: X|Y|Z #comment
    # if a "#" is found, then set 0 ! #WARNING: this is treated differently in Octopus. If norm(GaugeVectorField)==0, InducedFields is automatically enabled. 
    if ("#" in GaugeVectorField_tmp[0:5] or GaugeVectorField_inp2[0:2]): #NOTE: GaugeVectorField_tmp contains the amplitudes. GaugeVectorField_inp2 contains the keyword declaration. 
        GaugeFieldAmplitude = 0e0
    else: 
        # No "#" can be found in GaugeVectorField keyword, therefore Octopus switches to compute InducedFields without notifying. 
        try: 
            GaugeFieldAmplitude = GaugeVectorField_tmp[0] #returns the X occurrence
        except:
            GaugeFieldAmplitude = 0e0
            logger.warning((Header+"Warning: Selection of r.h.s. failed in INP. No data from "+str(GaugeVectorField_tmp)))
    if(GaugeFieldAmplitude==''):
        del GaugeFieldAmplitude
        GaugeFieldAmplitude = 0e0
    
    GaugeFieldDelay=0e0 #zero if was not used
    GaugeFieldDelay_tmp    = GaugeFieldDelay_inp.split("=")
    try:
        GaugeFieldDelay    = GaugeFieldDelay_tmp[1] #rhs
        GaugeFieldDelay    = float(GaugeFieldDelay) #if this fails, we set -1. 
    except: 
        GaugeFieldDelay=float('nan')
        logger.warning((Header+"Warning: Selection of r.h.s. failed in INP. No data in:"+str(GaugeFieldDelay_tmp)))
    
    ### Captures the functional which has been used for the computation.  #BUG: crashes when using TB09
    Functional_value = "lda" #default value: lda_x + lda_c_pz_mod
    try: 
      Functional_string = libString.LookUpStringInFile(filename, "XCFunctional")
    except:
      Functional_string = 'lda'
      logger.warning((Header+"** Warning: search of XCFunctional failed. Default value was set: lda"))
    
    if( len(Functional_string) > 0 ): #should be run only if Functional_string is defined. 
      Functional_value = Functional_string[13:] #capturing name of the functional
      Functional_lda_position = Functional_value.find("lda_x")
      Functional_gga_position  = Functional_value.find("GGA")
      Functional_mgga_position = Functional_value.find("MGGA")
      Functional_tb09_position = Functional_value.find("mgga_x_tb09")
      if(Functional_lda_position >= 0 and Functional_mgga_position == -1): #string was found
        Functional = IndexFunctional['LDA'] #lda
        logger.info((Header+"** Info: Functional LDA detected."))
      elif(Functional_mgga_position >= 0 and Functional_lda_position == -1): #MMGA was enabled found
        Functional = IndexFunctional['MGGA'] #mgga
        logger.info((Header+"** Info: Functional MGGA detected."))
      elif(Functional_gga_position >= 0 and Functional_mgga_position == -1 and Functional_lda_position == -1): #GGA was found, not MGGA.
        Functional = IndexFunctional['GGA'] #GGA
        logger.info((Header+"** Info: Functional GGA detected."))
      elif(Functional_lda_position == -1 and Functional_tb09_position >= 0):
        Functional = IndexFunctional['TB09'] #TB09
        logger.info((Header+"** Info: Functional TB09 detected."))
      else: 
        logger.error((Header+"** Error: XCFunctional contained unknown case."+str(Functional_value)))
        exit(-1)
    elif ( Functional_string == 'lda' ): #XCFunctional was not specified: assuming LDA.
      logger.warning((Header+"Warning: XCFunctional was not specified. Assuming default value: lda."))
      Functional = 0

    
    ### Captures if FreezeHXC was used. 
    FreezeHXC = 0
    try:
      FreezeHXC_string = libString.LookUpStringInFile(filename, "TDFreezeHXC")
      logger.info((Header+"Info: TDFreezeHXC_string: "+str(FreezeHXC_string)))
    except:
      logger.info((Header+"** Info: TDFreezeHXC was not found in the input file."))
      #FreezeHXC_string = ''
    if("yes" in FreezeHXC_string): #We shall not make a strict test, which was failing. 
      FreezeHXC = 1
      logger.info((Header+"Info: TDFreezeHXC enabled."))
    else:
      FreezeHXC = 0
      logger.info((Header+"Info: TDFreezeHXC not found."))
    
    # Captures if MoveIon was enabled. 
    MoveIon = 0
    try: 
      MoveIonString = libString.LookUpStringInFile(filename, "MoveIons")
      logger.info((Header+"Info: MoveIon_string: "+str(MoveIonString)))
    except:
      logger.info((Header+"Info: MoveIon was not found in the input file."))
    if("yes" in MoveIonString):
      MoveIon=1
    else:
      MoveIon=0
      logger.info((Header+"Info: MoveIon was not found in the input file."))
    
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

    logger.info((Header+"Wavelength (nm): "+str(wavelength_SI)))
    logger.info((Header+"Pulse duration (s): "+str(PulseDuration_SI)))
    logger.info((Header+"Double pulse delay (s): "+str(DoublePulseDelay_SI)))
    logger.info((Header+"Carrier Envelope Phase (deg): "+str(CEP1_SI * 180./pi)))
    
    # Catch the peak field #We could also take the one of input file..
    time                = np.asfarray(laser[:,1]) #column 2
    A1x                 = np.asfarray(laser[:,2]) #column number 3
    A2x                 = np.asfarray(laser[:,5]) #column number 6
    timestepTable       = np.gradient(time) #size N-1
    TimeStep            = timestepTable[0] #they are all equal, just select the first
    Efield              = - (np.gradient(A1x)/np.gradient(time)) #correct value is *1/c 
    Efield2             = - (np.gradient(A2x)/np.gradient(time))
    Efield_AU           = Efield/au.c_AU #fixing dimension #TODO: clean me that...
    Efield2_AU          = Efield2/au.c_AU
    Efield_amplitude_SI = np.max(abs(au.Field_AU_to_SI(Efield_AU)))
    Efield2_amplitude_SI= np.max(abs(au.Field_AU_to_SI(Efield2_AU)))
    
    logger.info((Header+"Electric field amplitude 1 = "+str(Efield_amplitude_SI / 1E9)+" V/nm."))
    logger.info((Header+"Electric field amplitude 2 = "+str(Efield2_amplitude_SI / 1E9)+" V/nm."))
    
    # Capturing cell volume
    try:
         SuperCellVolume_AU = calcVol.VolumeOfSupercell_Log(FolderName+"../")
    except:
         logger.warning((Header+"** Warning: sample computation volume was jumped!"))
         continue
    SuperCellVolume_SI = au.Volume_AU_to_SI(SuperCellVolume_AU)
    SampleCellVolume_AU = SuperCellVolume_AU #- calcVol.VolumeOfVacuum(FolderName+"../") #Volume of the sample (without vacuum region) # WARNING: we disabled the SLAB volume computation !
    CellVolume_AU = calcVol.VolumePerPoint_Log(FolderName+"../") #0.06385 #0.07820 #Bohr**3
    SampleCellVolume_SI = au.Volume_AU_to_SI(SampleCellVolume_AU)
    
    logger.info(Header+"** Sample volumes **")
    logger.info(Header+"Super-cell total volume = "+str(SuperCellVolume_AU)+" Bohr**3.")
    logger.info(Header+"Super-cell total volume = "+str(SuperCellVolume_SI)+" m^{3}.") 
    logger.info(Header+"")
    logger.info(Header+"Sample volume = "+str(SampleCellVolume_AU)+" Bohr**3.")
    logger.info(Header+"Sample vacuum volume = "+str(SampleCellVolume_SI)+" m^{3}.")
    logger.info(Header+"")
    logger.info(Header+"Numerical cell volume = "+str(CellVolume_AU)+" Bohr**3.")
    logger.info(Header+"Numerical cell volume = "+str(au.Volume_AU_to_SI(CellVolume_AU))+" m^{-3}.")
    logger.info(Header+"") 
    
    logger.info((Header+"** Importing simulation data..."))
    
    # Catch the final excited electron density. 
    try: 
        timeNe           = np.asfarray(Ne[:,1]) #WARNING: this does exist for metals...
    except: 
        logger.error((Header+"** Error: It seems this run was not complete. We jump it."))
        continue
    Nproj_AU         = np.asfarray(Ne[:,2]) #number of particles, not a density!
    density_proj_AU  = Nproj_AU / SampleCellVolume_AU #This is now a density (a.u.)
    density_proj_SI  = au.Density_AU_to_SI(density_proj_AU)
    sizeOfDensityTab = density_proj_SI.size

    density_proj_SI_final         = density_proj_SI[sizeOfDensityTab-1] ##NOTE: sometimes, n_ex changes with time after the pulse (particularly at low densities. This can change the slope of dN_exc/dI, inducing errors on necessary number of photons. This could explain the weird variations found in various results (e.g. Apostolova). It is better to capture the value right after the pulse. 
    density_proj_SI_max           = density_proj_SI.max()
    density_irreversibility_ratio = density_proj_SI_final / density_proj_SI_max
    
    logger.info((Header+"** Electron density analysis "))
    logger.info((Header+"Final electron density = "+str(density_proj_SI_final)))
    logger.info((Header+"Max.  electron density = "+str(density_proj_SI_max)))
#    print Header+"Irreversibility ratio (Final / Max) = "+str(density_irreversibility_ratio)
    ##==== Prepare plotting of energies ====
    timeEnergies   = np.asfarray(Energy[:,1]) #atomic units (! 1 is 2nd column)
    TotalEnergy    = np.asfarray(Energy[:,2]) #Hartree
    IonEnergy      = np.asfarray(Energy[:,3]) #Hartree
    ElectronEnergy = np.asfarray(Energy[:,5]) #Hartree
    Eigenvalues    = np.asfarray(Energy[:,6]) #Hartree
    Hartree        = np.asfarray(Energy[:,7]) #Hartree
    Exchange       = np.asfarray(Energy[:,9]) #Hartree
    
    ## NOTE: Octopus does not compute Hartree energy at all timesteps, it would be too expensive. 
    ## Redefining the zero with interpolated values between defined points. 
    Hartree_raw = Hartree
    Hartree_nonzero = np.nonzero(Hartree_raw) #returns index of elements that are non zero, but a strange dimension
    Hartree_nonzero = np.asarray(Hartree_nonzero)[0] #correcting dimension problem. flatnonzero() could also do it in principle.
    # print("** Debug: shape Hartree: "+str(np.shape(Hartree)))
    # print("**        shape Hartree_nonzero: "+str(np.shape(Hartree_nonzero)))
    Hartree_interp = interp1d(timeEnergies[Hartree_nonzero], Hartree[Hartree_nonzero], fill_value="extrapolate") #this remains to be a function
#    Hartree_interp = UnivariateSpline(timeEnergies[Hartree_nonzero], Hartree[Hartree_nonzero], k=1)
    # print("**        shape Hartree_interp: "+str(np.shape(Hartree_interp))) 
    Hartree        = Hartree_interp(timeEnergies)
    # except ValueError:
    #     print(Header+"Hartree energy could not be interpolated. I kept the initial value with ZERO when octopus was lazy to compute.")
    #     Hartree        = Hartree_nonzero #WARNING: not effective. 
#    del Hartree_interp, Hartree_nonzero
    
    ## NOTE: Octopus does not compute Exchange-Correlation energy at all timesteps, it would be too expensive. 
    ## Redefining the zero with interpolated values between defined points. 
    Exchange_nonzero = np.nonzero(Exchange)
    try: #tb09 keeps exchange == 0
        Exchange_interp  = interp1d(timeEnergies[Exchange_nonzero], Exchange[Exchange_nonzero])
        Exchange         = Exchange_interp(timeEnergies)
        del Exchange_interp, Exchange_nonzero
    except: 
        logger.warning((Header+"** Warning: could not interpolate energy output."))
    
    # Normalizing the database energies
    
    # We capture the energy increase due to the laser field. 
    # But there exist cases where the oscillations of the energy 
    # are larger than the energy amplitude. Therefore, lets just remove the ground state value. 
    ElectronEnergy = ElectronEnergy - ElectronEnergy[0]
    TotalEnergy    = TotalEnergy    - TotalEnergy[0]
    IonEnergy      = IonEnergy      - IonEnergy[0]
    Hartree        = Hartree        - Hartree[0]
    Exchange       = Exchange       - Exchange[0]
    Eigenvalues    = Eigenvalues    - Eigenvalues[0]
    
    # Extracting the FINAL shift current (!!!)
    timeCurrent           = np.asfarray(TotalCurrent[:,1]) #atomic units
    TotalCurrent_X_AU     = np.asfarray(TotalCurrent[:,2]) #atomic units
    TotalCurrent_Y_AU     = np.asfarray(TotalCurrent[:,3]) #a. u.
    TotalCurrent_Z_AU     = np.asfarray(TotalCurrent[:,4]) #a. u. 
    TotalCurrent_size     = TotalCurrent_X_AU.size
 
    timeCurrent_SI        = au.Time_AU_to_SI(timeCurrent)
    TotalCurrent_X_SI     = au.Current_AU_to_SI(TotalCurrent_X_AU)
    TotalCurrent_Y_SI     = au.Current_AU_to_SI(TotalCurrent_Y_AU) 
    TotalCurrent_Z_SI     = au.Current_AU_to_SI(TotalCurrent_Z_AU) 

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

    # TotalCurrent_SI_norm  = np.sqrt( np.power(TotalCurrent_X_SI,2)+np.power(TotalCurrent_Y_SI,2)+np.power(TotalCurrent_Z_SI,2) ) #computes norm of the current
 
    # Computing External field reconstructed using the current. Work in atomic units. 
    minLength=np.min([len(ElectronEnergy)-2, len(Efield_AU)-2]) ## BUGFIX: this is a nasty bugfix. It is highly unclear why dimension of Energy can be LONGER than the field! 

    # Efield_AU_shrinked=Efield_AU[:len(ElectronEnergy)-1] #will be used later

    Ein_AU     = Efield_AU[:minLength]
    PolarizationField_AU = cumtrapz(TotalCurrent_X_AU,x=timeCurrent)
    PolarizationField_AU_checkint = np.cumsum(TotalCurrent_X_AU)*TimeStep #integral of the current in time

    Ein_SI     = au.Field_AU_to_SI(Ein_AU)
    PolarizationField_SI = cumtrapz(TotalCurrent_X_SI[:minLength],x=timeCurrent_SI[:minLength])
    PolarizationField_SI_checkint = np.cumsum((TotalCurrent_X_SI[:minLength])*au.Time_AU_to_SI(TimeStep))

    Eext_AU    = Ein_AU + 4.*pi*PolarizationField_AU[:minLength] #builds E_ext(t) = E_in(t)+P(t)/epsilon_0 

    EextMinus_AU=Ein_AU - 4.*pi*PolarizationField_AU[:minLength] #builds E_ext(t) = E_in(t)+P(t)/epsilon_0

    Eext_SI    = Ein_SI[:minLength-1] + PolarizationField_SI[:minLength]/epsilon_0 
    EextMinus_SI= Ein_SI[:minLength-1] - PolarizationField_SI[:minLength]/epsilon_0

    Eext_AU_max=Eext_AU[minLength//2] #this gives the amplitude of the corresponding external field 
    EextMinus_AU_max=EextMinus_AU[minLength//2]
    Eext_SI_max=Eext_SI[minLength//2] #this gives the amplitude of the corresponding external field
    EextMinus_SI_max=EextMinus_SI[minLength//2] #this gives the amplitude of the corresponding external field

    logger.info((Header+"** Info: maximum external field +P (a.u.): "+str(float(Eext_AU_max))))
    logger.info((Header+"** Info: maximum external field -P (a.u.): "+str(float(EextMinus_AU_max))))
 
    logger.info((Header+"** Info: maximum external field +P (V/m): "+str(float(Eext_SI_max))))
    logger.info((Header+"** Info: maximum external field -P (V/m): "+str(float(EextMinus_SI_max))))
    
    # Computing integral(- J(t).E(t) dt)
    EnergyPoynting_AU = - np.multiply(Ein_AU[:minLength], TotalCurrent_X_AU[:minLength]) #J(t).E(t) in a.u.
    
    # Sept 2020: Lets get rid of the volume here
    #EnergyPoynting_AU_sum = SampleCellVolume_AU*cumtrapz(EnergyPoynting_AU[:minLength], x=timeCurrent[:minLength]) #int(J.E, dt) (Hartree) # WARNING: in plotTimeKeldysh.py, we didnt multiply by the volume. J(t) already contains the volume and field is homogeneous... This must be a mistake. 
    EnergyPoynting_AU_sum = cumtrapz(EnergyPoynting_AU[:minLength], x=timeCurrent[:minLength]) #int(J.E, dt) (Hartree)
    EnergyPoynting_AU_sum_checkint = np.cumsum(EnergyPoynting_AU[:minLength])*TimeStep #int(J.E, dt) (Hartree)

    EnergyPoynting_AU_sum_max = np.max(np.max(EnergyPoynting_AU_sum)) #maximum absorbed energy (in Hartree)
    EnergyPoynting_AU_sum_final = EnergyPoynting_AU_sum[len(EnergyPoynting_AU_sum)-1] #final absorbed energy #NOTE: false result!
 
    # Computing integral(J(t).E(t) dt)
    EnergyPoynting_SI = - np.multiply(Ein_SI[:minLength], TotalCurrent_X_SI[:minLength]) #J.E in a.u.
    # EnergyPoynting_SI_sum = SampleCellVolume_SI/e*cumtrapz(EnergyPoynting_SI[:minLength], x=timeCurrent_SI[:minLength]) #int(J.E, dt) [J/m3] * Vol[m3] = Joules. #BUG: Integration on too small numbers? WARNING Divergenc²    EnergyPoynting_SI_sum_max = np.max(np.max(EnergyPoynting_SI_sum)) #maximum absorbed energy
    EnergyPoynting_SI_sum = cumtrapz(EnergyPoynting_SI[:minLength], x=timeCurrent_SI[:minLength])/e #int(J.E, dt) #Volume was removed since it must be already contained in J(t).
    #EnergyPoynting_SI_sum_checkint = SampleCellVolume_SI/e*np.cumsum(EnergyPoynting_SI[:minLength])*au.Time_AU_to_SI(TimeStep)
    EnergyPoynting_SI_sum_checkint = np.cumsum(EnergyPoynting_SI[:minLength])*au.Time_AU_to_SI(TimeStep)/e #volume was removed here
    
    EnergyPoynting_SI_sum_final = EnergyPoynting_SI_sum[len(EnergyPoynting_SI_sum)-1] #final absorbed energy
   
    # EnergyPoynting_SI_check           =au.Energy_Hartree_to_eV(EnergyPoynting_AU)
    # EnergyPoynting_SI_sum_check       =au.Energy_Hartree_to_eV(EnergyPoynting_AU_sum)
    # EnergyPoynting_SI_sum_max_check   =au.Energy_Hartree_to_eV(EnergyPoynting_AU_sum_max)
    EnergyPoynting_SI_sum_final_check =au.Energy_Hartree_to_eV(EnergyPoynting_AU_sum_final)
#    print(Header+'** Info: max absorbed Poynting energy (Hartree): '+str(EnergyPoynting_AU_sum_max))
#    print(Header+'** Info: max absorbed Poynting energy (eV): '+str(EnergyPoynting_SI_sum_max))
#    print(Header+'** Info: max absorbed Poynting energy CHECK (eV): '+str(EnergyPoynting_SI_sum_max_check))
    logger.info("\n")
#    print(Header+"** Debug: timestep (AU): "+str(timeCurrent[1]-timeCurrent[0]))
#    print(Header+"** Debug: timestep (AU): "+str(TimeStep))
    logger.debug((Header+"** Debug integral (AU): "+str(np.max(np.max(PolarizationField_AU)))+", "+str(np.max(np.max(PolarizationField_AU_checkint))))) #VALID: well converges when sufficiently large
    logger.debug((Header+"** Debug integral (SI): "+str(np.max(np.max(PolarizationField_SI)))+", "+str(np.max(np.max(PolarizationField_SI_checkint))))) #NOTE: well converges when sufficiently large

    logger.debug((Header+"** Debug: Energy Poynting integration: "+str(np.max(np.max(EnergyPoynting_AU_sum)))+", "+str(np.max(np.max(EnergyPoynting_AU_sum_checkint))))) #NOTE: converges well when sufficiently large.
    logger.debug((Header+"** Debug: Energy Poynting integration (SI): "+str(np.max(np.max(EnergyPoynting_SI_sum)))+", "+str(np.max(np.max(EnergyPoynting_SI_sum_checkint)))))
    logger.debug("\n")

    logger.info((Header+'** Info: final absorbed Poynting energy (Hartree): '+str(EnergyPoynting_AU_sum_final)))
    logger.info((Header+'** Info: final absorbed Poynting energy (Hartree->eV): '+str(EnergyPoynting_SI_sum_final_check)))
    logger.info((Header+'** Info: final absorbed Poynting energy (eV): '+str(EnergyPoynting_SI_sum_final)))
    
    logger.info(Header+"Final DC current (SI): X: "+str(float(au.Current_AU_to_SI(TotalCurrent_X_AU_final)))+", Y: "+str(float(au.Current_AU_to_SI(TotalCurrent_Y_AU_final)))+", Z: "+str(float(au.Current_AU_to_SI(TotalCurrent_Z_AU_final)))+" C/m²/s.")
    TotalCurrent_SI_final = abs(TotalCurrent_X_AU_final) + abs(TotalCurrent_Y_AU_final) + abs(TotalCurrent_Z_AU_final)
    logger.info(Header+"** Final DC current: "+str(float(au.Current_AU_to_SI(TotalCurrent_SI_final)))+" (S.I.) ")
    
    # Extracting the final energy increase (total, and electronic)
    # NOTE: Total energy is constant with time. Then we use total energy field for Hartree term for now. In future, we may add a field for kinetic energy of atoms. 
    HartreeEnergy_AU_final = Hartree[Hartree.size-1] #WARNING: this may be equal to 0 most of the time. This could explain the discontinuity. 
    HartreeEnergy_SI_final = float(au.Energy_Hartree_to_eV(HartreeEnergy_AU_final))
    
    ExchangeEnergy_AU_final = Exchange[Exchange.size-1]
    ExchangeEnergy_SI_final = float(au.Energy_Hartree_to_eV(ExchangeEnergy_AU_final))
    
    EigenvaluesEnergy_AU_final = Eigenvalues[Eigenvalues.size-1]
    EigenvaluesEnergy_SI_final = float(au.Energy_Hartree_to_eV(EigenvaluesEnergy_AU_final))
    
    ElectronEnergy_AU_final = ElectronEnergy[ElectronEnergy.size-1]  # WARNING: we shall select the first moment after the laser pulse comes back to zero [because induced fields can actually decrease the energy AFTER the pulse, for obscure reasons.
    ElectronEnergy_SI_final = float(au.Energy_Hartree_to_eV(ElectronEnergy_AU_final))
    logger.info(Header+"End electron energy   increase (eV)     = "+str(ElectronEnergy_SI_final))
    logger.info(Header+"End Hartree  energy   increase (eV)     = "+str(HartreeEnergy_SI_final))
    logger.info(Header+"End exchange energy   increase (eV)     = "+str(ExchangeEnergy_SI_final))
    logger.info(Header+"End eigenvalue energy increase (eV)     = "+str(EigenvaluesEnergy_SI_final))
    logger.info(Header+"Final electron density = "+str(density_proj_SI_final))

    ## HERE WE TRY TO CAPTURE INSTANT AFTER THE PULSE. 
    # Several conditions should be simultaneously satisfied (classified by most restrictive conditions)
    # 2. t > tau_1+tau_2+delay, where tau_1,2 depend on the pulse shape. SinSquare: tau_1,2->2*tau_1,2. STH: tau_1 -> tau_1+2*2*pi/omega_1 (idem for tau_2). 
    # 1. E_tot^laser ~ 0 along with dE_tot^laser/dt~0. This also selects time before the laser pulse, and the intermediate region if two pulses are used. 
    # 3. At least a slow variation of n_e and energy: dn_e/dt~0 and dE_t/dt~0. 
    if(True):
#    try: 
        logger.info(Header+"** Info: trying condition for capturing ElectronDensity n_exc and ElectronEnergy right after the pulse: ")
        tolerance_field_nemesh=1E-3; tolerance_Efield = tolerance_field_nemesh #1E-3 #atomic units
        tolerance_Efield_relative=1E-7 #tolerance for selecting the end of bicolor pulses. Some of them have non varying density anymore, particularly during the last laser oscillations. 

        # Accounting for a second field, then only considering differential of field.
        Efield_12_AU=np.maximum(np.abs(Efield_AU),np.abs(Efield2_AU))
        # Computing differential of Efield_12_AU (then we can transparently work with two fields). 
        Efield_12_AU_diff=np.gradient(Efield_12_AU) #dE.

        # when using induced fields). Preparing dNe, dEnergy, d(Efield1)
        Efield_12_AU_max=np.amax(np.amax(np.abs(Efield_12_AU))) #maximum amplitude of the field (AU)
        FieldRelative = np.abs(Efield_12_AU)/Efield_12_AU_max #meaningless: we should compare distance to zero, normalized
        FieldRelativeVariation = np.abs(Efield_12_AU_diff)/Efield_12_AU_max #time-derivative dE_1+2, compared relatively to maximum field. dE/E should be small.  
        FieldRelativeVariation_bool_t = (FieldRelativeVariation < tolerance_Efield_relative)
#        if(not FieldRelativeVariation_bool_t[:len(FieldRelativeVariation_bool_t)]): 
        FieldEqualZero=(FieldRelative < tolerance_Efield_relative) #atomic units #this provides an array indexed on timestep where field does no longer variates (accounts for double pulse cases)
        FieldDiffEqualZero = FieldRelativeVariation_bool_t #(Efield_12_AU_diff < tolerance_Efield) 
        
        logger.debug((Header+"** Debug: Shape(FieldEqualZero    )="+str(np.shape(FieldEqualZero))))
        logger.debug((Header+"** Debug: Shape(FieldDiffEqualZero)="+str(np.shape(FieldDiffEqualZero))))

        # CONDITION 2. Time instants should be large enough. WARNING: octopus input data are in SI (our fault),
        # octopus files are in atomic units - we want to perform selection in atomic units to avoid roundoff).
        # Preparing tau_1, tau_2 definitions as function of pulse shape. 
        if(PulseTypeName=="SinSquared"):
            duration_tau_1 = 2*PulseDuration_SI  #SI
            duration_tau_2 = 2*PulseDuration2_SI #SI
        elif(PulseTypeName=="StH"):
            omega_1 = 2*pi*c/wavelength_SI  #SI
            omega_2 = 2*pi*c/wavelength2_SI #SI
            duration_tau_1 = PulseDuration_SI  + RiseTime*pi/omega_1 #SI
            duration_tau_2 = PulseDuration2_SI + RiseTime*pi/omega_2 #SI
        else: 
            logger.error(Header+"** ERROR: capture of end of pulse is not managed. Pulse type: "+str(PulseTypeName))
            raise IndexError #TODO: prepare Exception("NotImplemented")
            # raise IndexError #TODO: prepare Exception("NotImplemented")
        # Computing a time right after the pulse, accounting for situations of double pulse of SinSquared or StH shapes. 
        # The criterion of selection slightly varies if we use single or double pulses. It is 
        # more tolerant for double pulses (it is meant to re-use runs that are sufficiently converged in density, although the double pulse is not completely finished).   
        timeFarEnough = au.Time_SI_to_AU(duration_tau_1) #AU #VALID
        if(Efield2_SI_input > 1.): #as was well used line 250 (TODO: make a function)
            timeFarEnough += duration_tau_2 + au.Time_SI_to_AU(DoublePulseDelay_SI) #AU ScalarValue
            SizeOfWorkingArray=np.amin([len(FieldRelativeVariation_bool_t), len(timeEnergies)])
            Mask_energy_1 = FieldRelativeVariation_bool_t[:SizeOfWorkingArray] #for bicolor pulses, conditions are softer, as we have actually valid computations to show there - where n_exc was converged, but not E_1+E_2(t)
            logger.warning("Maximum dE/E (normalized): "+str(np.amax(np.amax(FieldRelativeVariation))))
            logger.warning("Last dE/E (normalized): "+str(FieldRelativeVariation[len(FieldRelativeVariation)-1]))
            logger.warning("Relative variation of field at last cell is: "+str(FieldRelativeVariation))
        else: #singlecolor case
            SizeOfWorkingArray=np.amin([len(FieldEqualZero), len(timeEnergies)]) # limiting size to smallest file
            
        ## MASK FOR JOINT CONDITION 1 and CONDITION 2. 
        # Build a mask to select region with all the conditions together, applicable on any object expressed on time grid of laser (and energy, provided they are the same)
        # BUILDING MASK FOR ENERGY - does not require interpolartion
            logger.debug((Header+"** Debug: size of FieldEqualZero, FieldDiffEqualZero, (timeEnergies > timeFarEnough): "+str(FieldEqualZero)+", "+str(FieldDiffEqualZero)+", "+str((timeEnergies > timeFarEnough))))
            Mask_energy_1= np.logical_and(FieldDiffEqualZero[:SizeOfWorkingArray], FieldEqualZero[:SizeOfWorkingArray]) #E==0 and dE==0. 
        
        
        # ElectronEnergy_diff = np.diff(ElectronEnergy)  #this provides variation of energy with timeLaser grid
        # Ne_diff             = np.diff(density_proj_SI) #this provides variation of density with timeNe grid
        Mask_energy = np.logical_and(Mask_energy_1, (timeEnergies[:SizeOfWorkingArray] > timeFarEnough)) #Mask for the energy observable expressed in the timeLaser grid (provided laser and enery time grids are the same)
        logger.debug((Header+"** Debug: size of Mask for Energy: "+str(len(Mask_energy))))
        logger.debug((Header+"** Debug: Mask = "+str(Mask_energy)))
 
        logger.debug((Header+"** Debug: size of Efield_AU, Efield2_AU: "+str(len(Efield_AU))+", "+str(len(Efield2_AU))))
        logger.debug((Header+"** Debug: shrinked size of Efield_AU, Efield2_AU, Efield_12_AU: timeNe: "+str(len(Efield_AU))+", "+str(len(Efield_12_AU))+", "+str(len(Efield2_AU))+", "+str(len(timeNe))))
        logger.debug(Header+"** Debug: size of dEfield_12_AU, Efield2_AU, Efield_12_AU: timeNe: "+str(len(Efield_AU))+", "+str(len(Efield_12_AU))+", "+str(len(Efield2_AU)))
        logger.debug(Header+"** Debug: size of shrinked dEfield_12_AU, Efield2_AU, Efield_12_AU: "+str(len(Efield_AU))+", "+str(len(Efield_12_AU))+", "+str(len(Efield2_AU)))
        logger.debug((Header+"** Debug: maximum existing instants for time, timeNe "+str(np.max(time))+", "+str(np.max(timeNe))))
        logger.debug((Header+"** Debug: timeNe ="+str(timeNe)))
        logger.debug((Header+"** Debug: time   ="+str(time)))
        logger.debug((Header+"** Debug: Efield_12_AU="+str(Efield_12_AU)))
        logger.debug((Header+"** Debug: max(Efield_12_AU)="+str(np.max(Efield_12_AU))))

        # BUILDING MASK FOR DENSITY - requires interpolation
        # One has to interpolate a field function on same temporal grid than Ne
        Efield_AU_interp_func = interp1d(time[:len(Efield_12_AU)], Efield_12_AU) #fails when max(time of Efield) is equal to max(time of Ne), ie, when simulation reaches the end... I suggest to remove error raise upon interpolation on final point. 
        # Computing dEfield_AU_NeGrid using dEfield_timeGrid and interpolating (ie changing grid). It should be more precise. 
        
        logger.debug((Header+"** Debug: before defining interp dEfield_AU_interp_func: "+str(len(time))+", "+str(len(Efield_12_AU_diff)))) 
        dEfield_AU_interp_func = interp1d(time[:len(Efield_12_AU_diff)], Efield_12_AU_diff) #, fill_value="extrapolate") #DEBUG
        
        Efield_AU_interp_Ne   =  Efield_AU_interp_func(timeNe) #timeNe and timeLaser do have the same values, but not the same discretization, neither dimension
        Efield_AU_interp_error=  (np.max(np.max(Efield_12_AU))-np.max(np.max(Efield_AU_interp_Ne)))/np.max(np.max(Efield_12_AU))
        logger.debug((Header+"** Debug: Efield_AU interpolation: "+str(np.max(np.max(Efield_AU_interp_Ne)))+", not interpolated: "+str(np.max(np.max(Efield_12_AU)))+", error rel.: "+str(Efield_AU_interp_error)))

        logger.debug((Header+"** Debug: timeSource: "+str(np.max(time[:len(ElectronEnergy)]))+", timeTarget:"+str(timeNe)))
        dEfield_AU_interp_Ne  = dEfield_AU_interp_func(timeNe)
        dEfield_AU_interp_error=(np.max(np.max(Efield_12_AU_diff))-np.max(np.max(dEfield_AU_interp_Ne)))/np.max(np.max(Efield_12_AU_diff))
        logger.debug((Header+"** Debug: d(Efield_AU): interpolated: "+str(np.max(np.max(dEfield_AU_interp_Ne)))+", not interpolated: "+str(np.max(np.max(Efield_12_AU_diff)))+", error rel.: "+str(dEfield_AU_interp_error)))

        # Checking that relative errors are less than 1 percent
        if(Efield_AU_interp_error > 0.05 or dEfield_AU_interp_error > 0.05):
            logger.warning((Header+"** WARNING: Interpolation error on field was greater than 5%!"))
#              continue

        # Identifying Efield == 0 and dEfield==0 on timeNe grid.
        ## NOTE: the selection works on a time filtering basis. We need to identify time instants that are relevant for the selection. But they must be expressed in the right time basis. Therefore, there is a lot of interpolation. 
        FieldEqualZero_Ne =(np.divide(Efield_AU_interp_Ne, np.max(Efield_AU_interp_Ne)) < tolerance_field_nemesh) #Selects when field is null on timeNe grid. Now based on bicolor field. #TODO: add a tolerance here ! Example shows 1E-4 would be fine. 
        dFieldEqualZero_Ne=(dEfield_AU_interp_Ne < tolerance_Efield )
        logger.debug((Header+"** DEBUG: shapes of Masks for Ne: "+str(np.shape(FieldEqualZero_Ne))+", "+str(np.shape(dFieldEqualZero_Ne))))

        # Building the mask that should apply to any observable based on time-grid of Ne
        ## CONDITION 3. Searching for region where energy (or density) stagnate (ie not a strict stabilization, e.g.,
        MaskNe1=np.logical_and(FieldEqualZero_Ne[:len(dFieldEqualZero_Ne)], dFieldEqualZero_Ne)
        MaskNe =np.logical_and(MaskNe1, (timeNe > timeFarEnough))

        ## CHECKING MASK INTEGRITY ... ie if manifolds contain at least 1 element (1 value of interest). 
        logger.debug((np.any(MaskNe)))
        if(not np.any(MaskNe)):
            logger.error((Header+"** ERROR: pulse may not be finished. Check n_ex!"))
            raise IndexError
        if(not np.any(Mask_energy)):
            logger.error((Header+"** ERROR: pulse may not be finished. Check energy file!"))
            raise IndexError
        
        ## THIS SEEMS FOR SOLVING DIMENSION PROBLEMS
#       print(Header, "** Info: timesteps where we found conditions to capture EnergyAfterPulse: ")
#        timeNe2      =timeNe[:len(timeNe)-2]
#       print(timeNe2)
        timeNe3      =timeNe[MaskNe]
        # print(timeNe3)
#         breakpoint()
        # Applying Mask on Energies
#        timeEnergies2=timeEnergies[:len(timeEnergies)-2]
        timeEnergies3=timeEnergies[Mask_energy] #filtering with condition given by Mask_energy
        
        # We defines masks for "sufficiently far instants after pulses" on the observables 
        FarEnoughInTime=(timeEnergies3 > timeFarEnough)
        FarEnoughInTime_Ne=(timeNe3 > timeFarEnough)

        #this gives the instant that should be used to filter signals
        # timeEnergies4=timeEnergies3[FarEnoughInTime]
        # timeNe4      =timeNe3[FarEnoughInTime_Ne]

        ##WARNING: if pulse is too short, manifold must be empty, this exception ValueError
        timeToFilterFrom       = timeFarEnough; #timeEnergies4[0]  rises.
        # timeToFilterFromDensity=timeNe4[0]

        # timeEnergiesCondition = ( timeEnergies > timeToFilterFrom        ) # this should define a mask that only keeps timesteps larger than laser pulse
        # timeDensityCondition  = ( timeNe       > timeToFilterFromDensity )

        # We filter ElectronEnergy using these two conditions: 
        ElectronEnergyRightAfterPulse = ElectronEnergy[Mask_energy][0] ## WARNING: it seems temporal filter was not applied !!!
        ElectronDensityRightAfterPulse= density_proj_SI[MaskNe][0] 
        logger.info((Header+"** Info: energy right after the pulse: "+str(ElectronEnergyRightAfterPulse)))
        logger.info((Header+"** Info: density right after the pulse: "+str(ElectronDensityRightAfterPulse)))           

    ## WE FETCH THE VALUE OF ENERGY RIGHT AFTER THE PULSE. 
        ElectronEnergy_SI_afterpulse   = float(au.Energy_Hartree_to_eV(ElectronEnergy[Mask_energy][0]))
        HartreeEnergy_SI_afterpulse    = float(au.Energy_Hartree_to_eV(Hartree    [Mask_energy][0]))
        ExchangeEnergy_SI_afterpulse   = float(au.Energy_Hartree_to_eV(Exchange   [Mask_energy][0]))
        EigenvaluesEnergy_SI_afterpulse= float(au.Energy_Hartree_to_eV(Eigenvalues[Mask_energy][0]))
    #except IndexError: #if some table was empty (e.g., simulation was not finished)
        #print(Header+"** ERROR: energy or density after pulse could not be captured.")
        #ElectronEnergy_SI_afterpulse   = 0e0
        #HartreeEnergy_SI_afterpulse    = 0e0
        #ExchangeEnergy_SI_afterpulse   = 0e0
        #EigenvaluesEnergy_SI_afterpulse= 0e0
        #ElectronDensityRightAfterPulse = 0e0

        #timeToFilterFrom = 0e0
    else: 
    #except IndexError: #ValueError: #interpolation creates some bugs when maximum time of Ne == maximum time of laser (ie, when a run finishes...)
        logger.error((Header+"** ERROR: energy or density after pulse could not be captured. Maybe pulse was not finished?"))
        ElectronEnergy_SI_afterpulse   = 0e0
        HartreeEnergy_SI_afterpulse    = 0e0
        ExchangeEnergy_SI_afterpulse   = 0e0
        EigenvaluesEnergy_SI_afterpulse= 0e0
        ElectronDensityRightAfterPulse = 0e0
        timeToFilterFrom = 0e0
    """
    except ValueError: #interpolation creates some bugs when maximum time of Ne == maximum time of laser (ie, when a run finishes...)
        logger.error((Header+"** ERROR: energy or density after pulse could not be captured for obscure reason. Maybe an error from interp1d() function at the boundary."))
        ElectronEnergy_SI_afterpulse   = 0e0
        HartreeEnergy_SI_afterpulse    = 0e0
        ExchangeEnergy_SI_afterpulse   = 0e0
        EigenvaluesEnergy_SI_afterpulse= 0e0
        ElectronDensityRightAfterPulse = 0e0
        timeToFilterFrom = 0e0
#        exit()
    """
    logger.info(Header+"Afterpulse electron energy   increase (eV)     = "+str(ElectronEnergy_SI_afterpulse))
    logger.info(Header+"Afterpulse Hartree  energy   increase (eV)     = "+str(HartreeEnergy_SI_afterpulse))
    logger.info(Header+"Afterpulse exchange energy   increase (eV)     = "+str(ExchangeEnergy_SI_afterpulse))
    logger.info(Header+"Afterpulse eigenvalue energy increase (eV)     = "+str(EigenvaluesEnergy_SI_afterpulse))
    logger.info(Header+"Afterpulse excited electron density (m-3)      = "+str(ElectronDensityRightAfterPulse))

#     EigenvaluesEnergy_SI_final
    # Append to the table
    # Structure: 
    # wavelength1, wavelength2, tau1, tau2, 
    # Epeak1, Epeak2, 
    # CEP1, CEP2, delay, 
    # kGrid, MeshSize, 
    # NeMax, NeEnd, density_irreversibility_ratio, EnergyMax
    # One line is one octopus run. O_o
     
#    print Header+"Checking integrity before adding to table..."
    #print float(runid), float(wavelength_SI), float(PulseDuration_SI)
    #print float(Efield1_SI_input), float(CEP1_SI), float(wavelength2_SI)
    #print float(PulseDuration2_SI), float(Efield2_SI_input), float(CEP2_SI)
    #print float(DoublePulseDelay_SI), float(CrossPolarization_angle), float(TheoryLevel)
    #print float(FieldInduced_enabled)
#    print "** Debug: GaugeFieldAmplitude: "+str(GaugeFieldAmplitude)
    #print float(GaugeFieldDelay)
    #print float(Functional),  float(FreezeHXC), float(MoveIon), float(density_proj_SI_max)
    #print float(density_proj_SI_final), float(density_irreversibility_ratio)
    #print float(TotalCurrent_SI_final), float(TotalCurrent_AU_max), float(TotalCurrent_AU_Xmax)
    #print float(TotalCurrent_AU_Ymax), float(TotalCurrent_AU_Zmax), float(FinalEnergy_SI_final)
    #print float(density_proj_SI_final/PulseDuration_SI), PulseType, RiseTime, float(timeNe[len(timeNe)-1])
      
    StackToAdd = [ float(runid), float(RunDate), kpointgrid, realmeshgrid, float(wavelength_SI), float(PulseDuration_SI),
              float(Efield1_SI_input), float(polarization_angle), float(CEP1_SI), float(wavelength2_SI),
              float(PulseDuration2_SI), float(Efield2_SI_input), float(CEP2_SI), 
              float(DoublePulseDelay_SI), float(CrossPolarization_angle), float(TheoryLevel), 
              float(FieldInduced_enabled), 
              float(GaugeFieldAmplitude), 
              float(GaugeFieldDelay), 
              float(Functional),  float(FreezeHXC), float(MoveIon), 
              float(density_proj_SI_max), float(density_proj_SI_final), float(density_irreversibility_ratio), 
              float(TotalCurrent_SI_final), float(TotalCurrent_AU_max), float(TotalCurrent_AU_Xmax), 
              float(TotalCurrent_AU_Ymax), float(TotalCurrent_AU_Zmax), float(ElectronEnergy_SI_final), float(EigenvaluesEnergy_SI_final), float(HartreeEnergy_SI_final), float(ExchangeEnergy_SI_final), 
              float(ElectronEnergy_SI_afterpulse), float(EigenvaluesEnergy_SI_afterpulse), float(HartreeEnergy_SI_afterpulse), 
              float(ExchangeEnergy_SI_afterpulse), float(EnergyPoynting_SI_sum_final_check), float(Eext_SI_max), float(ElectronDensityRightAfterPulse), 
              float(density_proj_SI_final/PulseDuration_SI), float(ElectronDensityRightAfterPulse/PulseDuration_SI),
              PulseType, RiseTime, float(timeNe[len(timeNe)-1]), 
              float(timeToFilterFrom) ]

    
    logger.info(Header+"Stack contains:")
    logger.info(StackToAdd)
    logger.info(Header+"** Info: Sizeof(PlottingTable)="+str(np.shape(PlottingTable))+".\nSizeof(StackToAdd)="+str(np.shape(StackToAdd)))
    
    logger.info((Header+"** INTEGRITY control: checking time derivative of latest frames"))

    # Take the latest frame of energy or n_exc, and check that we have finished the laser pulse. # WARNING: This assumes timestep for Ne = timestep for Energy = timeStep for Laser.
#    print(Header+"Debug: Length of Efield_AU: "+str(len(Efield_AU)))

#    try:
#    print(Header+"Debug: Length of ElectronEnergy: "+str(len(ElectronEnergy)))
    
    Efield_SI_final = au.Field_AU_to_SI(Efield_AU[minLength])
 #   except: 
  #      print(Header+"Debug: Length of density_proj_SI: "+str(len(density_proj_SI)))
   #     Efield_SI_final = au.Field_AU_to_SI(Efield_AU[len(density_proj_SI)-2])


    
   #     print(Header+"** Warning: unable to check if laser pulse was finished.")
   
    logger.info((Header+"   Efield at end of simulation V/m: "+str(Efield_SI_final)))
    ## NOTE: we check that derivative of energy or of n_ex is close to 0. If it is not, we jump the entry. 
    if(np.abs(Efield_SI_final) > 1.): 
        logger.warning((Header+" ** WARNING: THE RUN WAS NOT FINISHED. WE REMOVE IT FROM DB."))
        continue
    
    if(FromScratch): #We construct a new database from scratch
      try:
        PlottingTable = np.vstack(( PlottingTable, StackToAdd )) 
        logger.info(Header+"** Info: these data were added to the stack.")
      except:
        Header+"** Error: Failed to stack the database when looping on folders. Sizeof(PlottingTable)="+str(np.shape(PlottingTable))+". Sizeof(StackToAdd)="+str(np.shape(StackToAdd))
    else: #We only add the mising runs to the existing database
      # We have to check if this runs is not already available in already existing Database. 
      PreviousDB_ID = Database #Should be result from loadtxt. 
      #print PreviousDB_ID
      if int(runid) in PreviousDB_ID: #if this run was found in existing DB, we jump it. This formulation was successfully tested separately. 
          logger.info(Header+"** Info: ============================ THIS RUN WAS ALREADY EXISTING IN DB. Jumping it. =================")
          continue
      else: 
        PlottingTable = np.vstack(( PlottingTable, StackToAdd ))
        logger.info(Header+"** Info: these data were NEW and were ADDED to the stack.")
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
# TODO: this function should be replaced by Pandas usage.
def ExtractData(PlottingTable): #{{{
  runid                        = PlottingTable[:,0]
  RunDate                      = PlottingTable[:,1]
  kpointgrid                   = PlottingTable[:,2]
  realmeshgrid                 = PlottingTable[:,3]
  wavelengths1                 = PlottingTable[:,4]
  pulses1                      = PlottingTable[:,5]
  Efield1                      = PlottingTable[:,6]
  polarization_angle           = PlottingTable[:,7]
  CEP1                         = PlottingTable[:,8]
  wavelengths2                 = PlottingTable[:,9]
  pulses2                      = PlottingTable[:,10]
  Efield2                      = PlottingTable[:,11]
  CEP2                         = PlottingTable[:,12]
  Delay                        = PlottingTable[:,13]
  CrossPolarization            = PlottingTable[:,14]
  TheoryLevel                  = PlottingTable[:,15]
  FieldInduced                 = PlottingTable[:,16]
  GaugeFieldAmplitude          = PlottingTable[:,17]
  GaugeFieldDelay              = PlottingTable[:,18]
  Functional                   = PlottingTable[:,19]
  FreezeHXC                    = PlottingTable[:,20]
  MoveIon                      = PlottingTable[:,21]
  Ne_max                       = PlottingTable[:,22]
  Ne_final                     = PlottingTable[:,23]
  irreversibility              = PlottingTable[:,24]
  TotalCurrent_final           = PlottingTable[:,25]
  TotalCurrent_max             = PlottingTable[:,26]
  TotalCurrent_AU_Xmax         = PlottingTable[:,27]
  TotalCurrent_AU_Ymax         = PlottingTable[:,28]
  TotalCurrent_AU_Zmax         = PlottingTable[:,29]
  ElectronEnergy_max           = PlottingTable[:,30]
  EigenvaluesEnergy_max        = PlottingTable[:,31]
  HartreeEnergy_max            = PlottingTable[:,32]
  ExchangeEnergy_max           = PlottingTable[:,33]
  ElectronEnergy_afterpulse    = PlottingTable[:,34]
  EigenvaluesEnergy_afterpulse = PlottingTable[:,35]
  HartreeEnergy_afterpulse     = PlottingTable[:,36]
  ExchangeEnergy_afterpulse    = PlottingTable[:,37]
  EnergyPoynting_SI_sum_final  = PlottingTable[:,38]
  Eext_SI_max                  = PlottingTable[:,39]
  Ne_afterpulse                = PlottingTable[:,40]
  Wpi                          = PlottingTable[:,41]
  Wpi_afterpulse               = PlottingTable[:,42]
  PulseType                    = PlottingTable[:,43]
  RiseTime                     = PlottingTable[:,44]
  FinalTime                    = PlottingTable[:,45]
  timeToFilterFrom             = PlottingTable[:,46]

  return runid, RunDate, kpointgrid, realmeshgrid, wavelengths1, pulses1, Efield1, polarization_angle, \
      CEP1, wavelengths2, pulses2, Efield2, CEP2, \
      Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude,\
      GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, \
      Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, \
      TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, \
      ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, \
      EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, Wpi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom
#}}}

# TODO: this also can be replaced by Pandas usage
IndexTable = {'runid':                        0, 
              'rundate':                      1,
              'kpointgrid':                   2,
              'realmeshgrid':                 3,
              'wavelengths1':                 4,
              'pulses1':                      5,
              'Efield1':                      6,
              'polarization_angle':           7,
              'CEP1':                         8,
              'wavelengths2':                 9,
              'pulses2':                      10,
              'Efield2':                      11,
              'CEP2':                         12,
              'Delay':                        13,
              'CrossPolarization':            14,
              'TheoryLevel':                  15,
              'FieldInduced':                 16,
              'GaugeFieldAmplitude':          17,
              'GaugeFieldDelay':              18,
              'Functional':                   19,
              'FreezeHXC':                    20,
              'MoveIon':                      21,
              'Ne_max':                       22,
              'Ne_final':                     23,
              'irreversibility':              24,
              'Current_final':                25,
              'Current_max':                  26,
              'CurrentX_max':                 27,
              'CurrentY_max':                 28,
              'CurrentZ_max':                 29,
              'ElectronEnergy_final':         30,
              'EigenvaluesEnergy_max':        31,
              'HartreeEnergy_max':            32,
              'ExchangeEnergy_max':           33,
              'ElectronEnergy_afterpulse':    34,
              'EigenvaluesEnergy_afterpulse': 35,
              'HartreeEnergy_afterpulse':     36,
              'ExchangeEnergy_afterpulse':    37,
              'EnergyPoynting_final':         38,
              'Eext_SI':                      39,
              'Ne_afterpulse':                40,
              'Wpi':                          41,
              'Wpi_afterpulse':               42,
              'PulseType':                    43,
              'RiseTime':                     44,
              'FinalTime':                    45,
              'FilterTime':                   46
              }

## Counts the number of available laser parameters in a set of octopus results
# @param Database: array of octopus results. It can be filtered before using FilterDatabase functions. See libDatabase.py. 
def CountTheNumberOfLaserParameters(Database):#{{{
  runid, RunDate, kpointgrid, realmeshgrid, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, \
  Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, \
  MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, \
  TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, \
  ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, \
  ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, Wpi, Wpi_afterpulse, \
  PulseType, RiseTime, FinalTime, FilterTime = ExtractData(Database)
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


  logger.info("====== PULSE 1 =====")
  logger.info("** Efield1*           size = "+str(len(Efield1_uniq)))
  logger.info("** wavelength1*       size = "+str(len(wavelengths1_uniq)))
  logger.info("** pulse duration*    size = "+str(len(pulses1_uniq)))
  logger.info("** CEP1               size = "+str(len(CEP1_uniq)))

  logger.info("====== PULSE 2 =====")
  logger.info("** Efield2*           size = "+str(len(Efield2_uniq)))
  logger.info("** wavelength2*       size = "+str(len(wavelengths2_uniq)))
  logger.info("** pulse duration*    size = "+str(len(pulses2_uniq)))
  logger.info("** CEP2               size = "+str(len(CEP2_uniq)))

  logger.info("=== THEORY LEVEL ===")
  logger.info("** Number of theory levels = "+str(len(TheoryLevel_uniq)))
  logger.info("** Induced Fields          = "+str(len(FieldInduced_uniq)))
  logger.info("** Functionals             = "+str(len(Functional_uniq)))
  logger.info("** FreezeHXC               = "+str(len(FreezeHXC_uniq)))
  logger.info("** MoveIon                 = "+str(len(MoveIon_uniq)))
  
  logger.info("=== COUPLING EFFECTS ===")
  logger.info("** Delays             size = "+str(len(Delay_uniq)))
  logger.info("** CrossPolarization  size = "+str(len(CrossPolarization_uniq)))

  logger.info("====== RESULTS =====")
  logger.info("** Ne_final           size = "+str(len(Ne_final)))
  logger.info("** Irreversibility    size = "+str(len(irreversibility)))
  logger.info("** Total energy       size = "+str(len(ElectronEnergy_max)))
  
  return(wavelengths1_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses1_uniq, wavelengths2_uniq, CEP2_uniq,
         Efield2_uniq, pulses2_uniq)
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
        PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
        PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, Efield1*0.99, IndexTable['Efield1'])
        PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, Efield1*1.01, IndexTable['Efield1'])
        PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
        PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0., IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0.
        PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, CEP, IndexTable['CEP1']) #removes fields where CEP != 0. #NOTE: Can be done for Si-LDA.
        PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0., IndexTable['FreezeHXC']) #removes fields where FreezeHXC != 0.
        PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0., IndexTable['FieldInduced']) #removes fields where FreezeHXC != 0.
        PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 1., IndexTable['TheoryLevel']) #removes runs where TheoryLevel != dft.
        PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, MoveIon, IndexTable['MoveIon'])
        PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, PulseType, IndexTable['PulseType']) #keeps only pulses of type StH (3)
        PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, pulseDuration, IndexTable['pulses1']) #keeps only pulses of selected duration in argument
        if(RiseTime >= 1): #select by RiseTime only if it is defined
                PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, RiseTime, IndexTable['RiseTime'])
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
        PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works 
        PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, Efield1*0.99, IndexTable['Efield1'])
        PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, Efield1*1.01, IndexTable['Efield1'])
        PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
        PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, wavelength1, IndexTable['wavelengths1'])
#        PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0., IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0. #TODO: management of cross polarization is PROBLEMATIC. See libPlotting.py:213 CrossPolarization_Angle issue.

#        PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, 0., IndexTable['CEP1']) #removes fields where CEP < 0. Just for wrong buggy reasons.
        PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0., IndexTable['FreezeHXC']) #removes fields where FreezeHXC != 0.
        PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, InducedFields, IndexTable['FieldInduced'])
        PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 1., IndexTable['TheoryLevel']) #removes runs where TheoryLevel != dft.
        PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, MoveIon, IndexTable['MoveIon'])
        PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, PulseType, IndexTable['PulseType']) #keeps only pulses of type StH (3)
#        PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, pulseDuration, IndexTable['pulses1']) #keeps only pulses of selected duration in argument # NOT EFFICIENT
        PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, pulseDuration*1.1, IndexTable['pulses1'])
        PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, pulseDuration*0.9, IndexTable['pulses1'])
        if(RiseTime >= 1): #select by RiseTime only if it is defined
                PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, RiseTime, IndexTable['RiseTime'])
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
PulseTypeID     = {'gaussian': 0, 'SinSquared': 1, 'STH': 2, 'StH': 3, 'StH4': 4, 'SinFourth': 5}
PulseTypeID_rev = {0: 'gaussian', 1: 'SinSquared', 2: 'STH', 3: 'StH', 4: 'StH4', 5: 'SinFourth'} #reversed

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
    timeCurrent_SI = timeCurrent #* au.Time_AU_to_SI(1.)
    
    #FinalDampingInstant = tmax #timeCurrent_SI[maxSize/2] #this one might be more reasonable
    logger.info(Header+"** Info: HHG dampingTime: "+str(round(DampingTime*1E15))+" fs.")
    logger.info(Header+"** Info: HHG final instant: "+str(round(FinalDampingInstant*1E15))+" fs.")
    logger.info(Header+"** Info: omega = "+str(round(omega)))

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
  logger.info(len(BrutX))
  logger.info(len(BrutY))
  logger.info(len(BrutZ))
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
      #cmap = cm.get_cmap(name='seismic', lut=None)
      #print("Colormap:", cmap)
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
      
      ## Create a diverging colormap centered on zero.
      # d e f DefineGradientColorMap(Z):
        #print(Header+"** Debug:"+np.min(Z))
        #cvals  = [np.min(Z), 0, np.max(Z)]
        #colors = ["red","white","blue"]
   
        #norm=plt.Normalize(min(cvals),max(cvals))
        #tuples = list(zip(map(norm,cvals), colors))
        #cmap = mp.colors.LinearSegmentedColormap.from_list("", tuples)
      #return cmap
      
      plt.savefig(filename+'.eps')
      plt.savefig(filename+'.png')
      logger.info(Header+"** Info: saved "+filename+"eps|png.")
      if(ShowFinalPlot==1):
        plt.show()
      
      if(PlotGradient): #{{{
        logger.info(Header+"** Info: Plotting the gradient...")
        logger.debug(Header+"** Debug: Shape of I1: "+str(np.shape(E1_t)))
        logger.debug(Header+"** Debug: Shape of I2: "+str(np.shape(E2_t)))
        logger.debug(Header+"** Debug: Shape of the grad(interpolated field X): "+str(np.shape(Nexc_t_grad[0])))
        logger.debug(Header+"** Debug: Shape of the grad(interpolated field Y): "+str(np.shape(Nexc_t_grad[1])))
        
        fig2=plt.figure()
        plt.title(r"$\nabla_{E_1}$"+QuantityTitle)
        #plt.xlabel(r'$I_1$ (W/m$^{2}$)')
        #plt.ylabel(r'$I_2$ (W/m$^{2}$)')
        plt.xlabel(r'$E_1$ (V/m)')
        plt.ylabel(r'$E_2$ (V/m)')
        data_Z=Nexc_t_grad[0]
        #print(Header+"** Debug"+str(np.shape(data_Z)))
        #cmap=DefineGradientColorMap(data_Z)
        CS = plt.contourf(E1_t, E2_t, data_Z, cmap=cmap)
        plt.triplot(triang, lw=0.1, color='grey')
        plt.colorbar(CS)
        plt.tight_layout()
        plt.savefig(filename+'-GradientE1.eps')
        plt.savefig(filename+'-GradientE1.png')
        if(ShowFinalPlot==1):
            plt.show()
            
        fig3=plt.figure()
        plt.title(r"$\nabla_{E_2}$"+QuantityTitle)
        plt.xlabel(r'$E_1$ (V/m)')
        plt.ylabel(r'$E_2$ (V/m)')
        data_Z=Nexc_t_grad[1]
        #cmap=DefineGradientColorMap(data_Z)
        #CS = plt.contourf(E1_t, E2_t, data_Z, cmap=cmap)
        CS = plt.contourf(E1_t, E2_t, data_Z, cmap=cmap)
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
def plot2dHeatMap_HHG(BrutX, BrutY, BrutZ, QuantityTitle, filename, ShowFinalPlot=True, PlotGradient=False, PlotTriangles=True, xlabel='', ylabel='', cmap_ext='RdBu_r', number_of_levels=10, xlimit=[0,10], ylimit=[0,10], Egap=2.56*e, resolutionX=100, resolutionY=100): #{{{
#    cmap='Blues'
#    cmap='RdBu_r'
#    afmhot_r, hot, hot_r, bone_r, 
#    colormap="bone_r"

  logger.info(len(BrutX))
  logger.info(len(BrutY))
  logger.info(len(BrutZ))
  X=list(set(BrutX))
  Y=list(BrutY)

#  YY, XX = np.meshgrid(X, Y)
  # If the size is big enough, let's make a plot.
  #if(len(DensityMax2D) >= 2): #{{{
  if(len(BrutZ) >= 2): #{{{
      ## CONCEPT: initial mesh is plotted using a non-structured grid. 
      ##          while the gradient is computed after interpolation on
      ##          a structured mesh
      logger.info("Building the unstructured mesh triangulation...")
      triang = tri.Triangulation(BrutX, BrutY) #computes Delaunay triangulation
#      triang.set_mask(np.hypot(BrutX[triang.triangles].mean(axis=1),
#                               BrutY[triang.triangles].mean(axis=1))
#                      < 0.15)
      # We can project the triangular mesh on a orthogonal by performing interpolation as done here. https://matplotlib.org/3.1.0/gallery/images_contours_and_fields/triinterp_demo.html

      logger.info("Preparing refinement objects...")
      refiner = tri.UniformTriRefiner(triang)
      tri_refi, z_test_refi = refiner.refine_field(BrutZ, triinterpolator=None, subdiv=0) #For strange reasons, refining induces some very bad result.

#      logger.info("Preparing TriCubic interpolation...")
#      interp_cubic = tri.CubicTriInterpolator(tri_refi, z_test_refi, kind='geom') #min_E, geom

#      logger.info("Refining the field...")
#      tri_refi, z_test_refi = refiner.refine_field(BrutZ, interp_cubic, subdiv=3)

      logger.info("Interpolation on a structured grid to compute the gradient...")
     # E1_t, E2_t = np.meshgrid( np.linspace(np.min(BrutX), np.max(BrutX), resolutionX), np.linspace(np.min(BrutY), np.max(BrutY), resolutionY) ) #this one returns weird pics for HHG with CEP
      E1_t, E2_t = np.meshgrid( np.linspace(np.min(BrutX), xlimit, resolutionX), np.linspace(np.min(BrutY), ylimit, resolutionY) ) #let's try to limit the axis of energies (Y-axis)
      interp_cubic = tri.CubicTriInterpolator(tri_refi, z_test_refi, kind='geom') #min_E, geom
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
      for n in [1,2,3,4,5,6]:
        plt.axvline(x=Egap/e/n, linestyle=":")
        plt.axhline(y=1E6*h*c/(Egap/e/n), linestyle=':') 
        plt.axvline(x=n*Egap/e, linestyle=":")
        # plt.axhline(y=1E6*h*c/(Egap*n), linestyle=':') #not relevant

      # plt.savefig(filename+'.eps')
      plt.savefig(filename+'.png')
      logger.info(Header+"** Info: saved "+filename+".png.")
      if(ShowFinalPlot==1):
        plt.show()
      
      if(PlotGradient):
        logger.info(Header+"** Info: Plotting the gradient...")
        logger.debug(Header+"** Debug: Shape of I1: "+str(np.shape(E1_t)))
        logger.debug(Header+"** Debug: Shape of I2: "+str(np.shape(E2_t)))
        logger.debug(Header+"** Debug: Shape of the grad(interpolated field X): "+str(np.shape(Nexc_t_grad[0])))
        logger.debug(Header+"** Debug: Shape of the grad(interpolated field Y): "+str(np.shape(Nexc_t_grad[1])))
        
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
