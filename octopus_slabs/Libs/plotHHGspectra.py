#!/usr/bin/env python
#-*- coding: utf-8 -*-

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

## @package plotHHGspectra
# Collects the generated time-dependent octopus simulations spectra HHG.csv from all subfolders (can last minutes)
# Then prepares a picture using the collected data for each identified field amplitude

import glob
import os

from spp_extended_theory.Keldysh import libKeldysh

from spp_extended_theory.Libs import libDatabase
# from spp_extended_theory.Keldysh.libKeldysh import *

from scipy.interpolate import interp1d
from scipy.constants import e, h, hbar, c
from octopus_slabs.Libs import libPlotting

from matplotlib import rc
import matplotlib as mp
import numpy as np

import pandas as pd

from   octopus_slabs.Libs.libLogging         import init_logger

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':'14'})
## for Palatino and other serif fonts use:
#rc('font', **{'family':'serif', 'serif':['Palatino'], 'size':'16'})
rc('text', usetex=False) #True: Does not work on Draco. Just put False, then, and it will be smooth. 
mp.rcParams['legend.numpoints'] = 1

Header="[plotHHGspectra] "
logger = init_logger(__name__, verbose=True) #"plotFinalQuantities")

ShowFinalPlot              = 1
HHG_samples                = 6534 #6534: good!! 600: insufficient!!

path = "*/td.general/"
DatabaseFile="Database.csv"

# NOTE: two ways to capture the corresponding data can be employed. 
# 1. Looking the individual input file.
# 2. Or accessing the database using the ID of the run. Must be much more clever way.
# METHOD 1: accessing each subfolder. 
def ParseAllFolders(): #{{{
	for FolderName in glob.glob(path): #{{{ #NOTE: this takes everything available: not sorted.
	    print(Header+"** Info: entering folder: "+str(FolderName)+"\n")
	  
	    ### METHOD 1: rebuilding database without choosing the order of the sequence.  
	    filename = FolderName+"../inp"
	    Efield1_SI = libPlotting.GetParameterValueFromInputFile(filename, "EfieldAsi = ")
	    Efield2_SI = libPlotting.GetParameterValueFromInputFile(filename, "EfieldBsi = ")
	    PulseDuration_SI = libPlotting.GetParameterValueFromInputFile(filename, "tauAsi = ")
	    PulseDuration2_SI= libPlotting.GetParameterValueFromInputFile(filename, "tauBsi = ")
	    wavelength1_SI   = libPlotting.GetParameterValueFromInputFile(filename, "lambdaA = ")
	    wavelength2_SI   = libPlotting.GetParameterValueFromInputFile(filename, "lambdaB = ") 
   
	    print("Pulse duration 1 (s): "+str(PulseDuration_SI))
	    print("Field amplitude 2 (V/m): "+str(Efield2_SI))
	
	    # NOTE: all of these data were already collected in the Database. Fetch them from database.  
	    # If the run considers only single pulses
	    if(Efield2_SI > 1): #we only select runs that are with single pulses. 
	        # NOTE: We need multiple wavelengths, but one spectrum per wavelength.       
	
	        # List all the files HHG-spectrums.csv in subfolders
	        # try: 
	        HHG=np.loadtxt(FolderName+"../td.general/HHG-Spectrum.csv", dtype='str', skiprows=0)
	        # except:
	        #  print "** Error: Loading of HHG-Spectrum.csv failed." 
	        # Size of the spectrum? 
	        HHG=np.asfarray(HHG)
	        print("Size of HHG file: "+str(len(HHG)))
	#        print HHG
	#        print "Extracting 1 column: "
	#        print HHG[:,1]
	        # Interpolate it to the size of the others
	       
	        # Add it as a line in table
	        X=np.append(HHG[:,0], X) #HHG frequencies for 1 run
	        Z=np.append(HHG[:,1], Z) #HHG amplitudes  for 1 run
	        print(counter)
	        Y=np.append(counter*np.ones(len(HHG)), [Y]) #Y should contain the parameter we variate. If taking it all, we use a counter to represent them all as a simple stack for now.
	        print(Y)
	        # Prepare the final table
	        # Plot the intermediate picture
	        # plt.figure()
	        # plt.plot(X, Z)
	        # plt.show() 
	        # Add analytical Stark effect there
	        # Cleaning
        	del HHG
	        counter+=1
	    else: 
	        print("** Info: Double pulses are ignored.")
	#}}} #glob loop
#}}}
    

Database = np.loadtxt(DatabaseFile, skiprows=3)
## ======== METHOD 2: construct the sub-database of interest

def PlotHHGsinglePulsewithWavelength(Efield, tau=20E-15, CEP=0, MoveIon=0, PulseType=3, RiseCycles=2, LogScale=False): #{{{
    
    # Declaring variables
    X=[]; Y=[]; Z=[]
   
    # FilterSinglePulseDatabaseWithParameters_wavelengths(Database, Efield1, pulseDuration, CEP, PulseType=3, RiseTime=2)
    #Efield =7.5E8 #V/m
    #tau    =20e-15
    #CEP    =0
    #MoveIon=0
    #PulseType=3; RiseCycles=2
    Filtered1= libPlotting.FilterSinglePulseDatabaseWithParameters_wavelengths(Database, Efield, tau, CEP, PulseType, RiseCycles)

    # SORTING
    # One must sort here the DB as function of the sequence we would like to plot.
    # Sorted1=np.sort(Filtered1, axis=0)
    Sorted2=Filtered1[np.lexsort((Filtered1[:, libPlotting.IndexTable['Efield1']], Filtered1[:, libPlotting.IndexTable['wavelengths1']]))]
    # NOTE: how to manage the duplicated computations? We would like to keep longest computations. But may not remove duplicates. For now: manually disabling duplicated entries in DB. 
    # One could use np.unique() function, but we may loose the JobID then? 
    # Sorted3=libDatabase.FilterDatabaseRemove( #lets remove the shortest duplicated run

    # Extracting the data from filtered database
    try: 
        runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, Energy_max, Wpi, PulseType, RiseTime, FinalTime  = libPlotting.ExtractData(Sorted2)
    except:
        print("Extraction of sorted table failed. Maybe some elements are duplicated? Tip: keep the longest runs. Length(Sorted2)="+str(len(Sorted2)))
        exit()

    # DuplicatedEntries=(list(set(runid))==len(Filtered1)) #test if we have duplicated entries
    # if(DuplicatedEntries): #we keep only the longest runs
    #  print "There exists DUPLICATED ENTRIES (several iterations of same parameters. Please clean the database.csv flie in that respect, then restart."
    #  print DuplicatedEntries
    #  exit()  

    for JobID in np.arange(0,len(runid)):
        runid_int = int(runid[JobID])
        FileLookUp="octopus-td."+str(runid_int)+".log"
        ListFound=glob.glob("*/*/"+FileLookUp)
        runid_folder=os.path.split(os.path.abspath(ListFound[0]))[0] #assumes only 1 finding per loop
        print(runid_folder)
        try:
            HHG=np.loadtxt(runid_folder+"/td.general/HHG-Spectrum.csv", dtype='str', skiprows=0)
        except:
            print("** ERROR: Loading of "+runid_folder+"/td.general/HHG-Spectrum.csv failed.") 
            continue #just jump this run. Keep plotting the rest. 
        # Size of the spectrum? 
        HHG=np.asfarray(HHG)
        print("Size of HHG file: "+str(len(HHG)))

        if(len(HHG[:,0])!=HHG_samples): #NOTE: 600 is just the most employed number of samples for FFT in our computations. The size of HHG was like 600 all the time. Let's interpolate to that size if different.
            print("** Info: interpolating...")
#            HHG_func=InterpolatedUnivariateSpline(HHG[:,0], HHG[:,1]) #does not work very well. Sure, they are not many points in input. 
            HHG_func=interp1d(HHG[:,0], HHG[:,1], kind='cubic')
            HHG_omega_new=np.linspace(np.min(HHG[:,0]), np.max(HHG[:,0]), HHG_samples)
            Z_interp=HHG_func(HHG_omega_new)
            X=np.append(HHG_omega_new, X) #HHG frequencies for 1 run
            Z=np.append(Z_interp,      Z) #HHG amplitudes  for 1 run
            #          print counter
            Y=np.append(np.multiply(1E6*wavelengths1[JobID],np.ones(len(HHG_omega_new))), [Y]) #we are sure that JobID is increasing monotonically. 
#            Y=np.append(np.multiply(CEP1[JobID]*180./np.pi,np.ones(len(HHG_omega_new[:,0]))),[Y])
        else: #dont interpolate
            X=np.append(HHG[:,0],X)
            Z=np.append(HHG[:,1],Z)
            Y=np.append(np.multiply(1E6*wavelengths1[JobID],np.ones(len(HHG[:,0]))),[Y]) #TODO: we can use something else to plot, such as CEP or wavelengths2 !
#            Y=np.append(np.multiply(CEP1[JobID]*180./np.pi,np.ones(len(HHG[:,0]))),[Y])
        print(str(int(1E9*wavelengths1[JobID]))+" nm")
        print(str(h*c/wavelengths1[JobID]/e)+" eV")
        del HHG

    #### ======== COMMON PART 
    # Converts Y to array of reals
    Y = np.asfarray(Y)

    # Returns log10 of Z, but clean.
    Z_clipped=np.clip(Z, 1, np.max(Z))
    if(LogScale): 
        Z_log=np.log10(Z_clipped)
    # Z_contourf=np.reshape(Z_log, len(counter, len(Z_log)/counter))
        Z=Z_log #Z_clipped
    else:
        Z=Z_clipped

    xlabel=r"$\hbar \omega_{HHG}$ (eV)"
    ylabel=r"Wavelength $\lambda$ ($\mu$m)"
    # ylabel=r"Photon energy $\hbar \omega$ (eV)"
    title=r"$HHG(\omega_{probe}$, $\lambda_{pump})$"
    #    cmap='Blues'
#   cmap='RdBu_r' #null values generate serious noise in the picture
#   NICE: 'bone_r', 'afmhot_r'
#   TO TEST: ocean.
    # Not amazing: hot_r
    colormap='ocean_r'

    ShowFinalPlot=False
    PlotTriangles=False
    PlotGradient =False

    xlimit=[0,10]   #eV
    ylimit=[0,2.5]  #um
    numberoflevels=30
    Egap=2.56*e
#plot2dHeatMap(BrutX, BrutY, BrutZ, QuantityTitle, filename, ShowFinalPlot=True, PlotGradient=False, PlotTriangles=True, xlabel='', ylabel='', cmap_ext='afmhot_r')
    filename="HHGspectrum-"+str(int(Efield))+"-HHGamples"+str(HHG_samples)
    if(LogScale):
        filename+="LogScale"

    ## EXPORT to external file.
    print("** Info: Exporting to file "+filename+".csv")
    libDatabase.ExportToTxt(np.transpose([X, Y, Z]), filename+".csv")

    print("** Preparing plot...")
    libPlotting.plot2dHeatMap_HHG(X, Y, Z, title, filename, ShowFinalPlot, PlotGradient, PlotTriangles, xlabel, ylabel, colormap, numberoflevels, xlimit, ylimit, Egap)
    return 0
#}}}

## Prepares color plot HHG(omega_HHG, CEP) from database
# LogScale does not make sense here (we look at fluctuations)
def PlotHHGsinglePulsewithCEP(wavelength=800e-9, Efield=1E9, tau=20E-15, MoveIon=0, PulseType=3, RiseCycles=2,
                              InducedFields=0): #{{{
    global Database 
    # Declaring variables
    X=[]; Y=[]; Z=[]
   
    # FilterSinglePulseDatabaseWithParameters_wavelengths(Database, Efield1, pulseDuration, CEP, PulseType=3, RiseTime=2)
    #Efield =7.5E8 #V/m
    #tau    =20e-15
    #CEP    =0
    #MoveIon=0
    #PulseType=3; RiseCycles=2
    # CHECKING SIZE OF DATABASE
    libPlotting.CountTheNumberOfLaserParameters(Database)

    print("INFO: FILTERING...")
    Filtered1= libPlotting.FilterSinglePulseDatabaseWithParameters_CEP(Database, wavelength, Efield, tau, PulseType, 0, RiseCycles, InducedFields)
    Filtered1=np.asfarray(Filtered1) #lets use only numbers there.

    # CHECKING size of returned table
    print("INFO: size of returned table: "+str(np.size(Filtered1)))
    
    # SORTING
    # One must sort here the DB as function of the sequence we would like to plot.
    # Sorted1=np.sort(Filtered1, axis=0)
    Sorted2=Filtered1[np.lexsort((Filtered1[:, libPlotting.IndexTable['Efield1']], Filtered1[:,libPlotting.IndexTable['CEP1']]))]
    # NOTE: how to manage the duplicated computations? We would like to keep longest computations. But may not remove duplicates. For now: manually disabling duplicated entries in DB. 
    # One could use np.unique() function, but we may loose the JobID then? 
    # Sorted3=libDatabase.FilterDatabaseRemove( #lets remove the shortest duplicated run

    # Extracting the data from filtered database
    # try: 
#    runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, Energy_max, Wpi, PulseType, RiseTime, FinalTime  = libPlotting.ExtractData(Sorted2)
    runid, RunDate, kpointgrid, realmeshgrid, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, \
    pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, \
    Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, \
    TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, \
    HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, \
    HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse,  EnergyPoynting_SI_sum_final, Eext_SI_max, \
    Ne_afterpulse, Wpi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom = libPlotting.ExtractData(Sorted2)
    # except:
    #    print("Extraction of sorted table failed. Maybe some elements are duplicated? Tip: keep the longest runs. Length(Sorted2)="+str(len(Sorted2)))
    #    exit()

    # DuplicatedEntries=(list(set(runid))==len(Filtered1)) #test if we have duplicated entries
    # if(DuplicatedEntries): #we keep only the longest runs
    #  print "There exists DUPLICATED ENTRIES (several iterations of same parameters. Please clean the database.csv flie in that respect, then restart."
    #  print DuplicatedEntries
    #  exit()  

    for JobID in np.arange(0,len(runid)):
        runid_int = int(runid[JobID])
        FileLookUp="octopus-td."+str(runid_int)+"*.log"
        ListFound=glob.glob("*/*/"+FileLookUp)
        print(FileLookUp)
        runid_folder=os.path.split(os.path.abspath(ListFound[0]))[0] #assumes only 1 finding per loop
        print(runid_folder)
        try:
            HHG=np.loadtxt(runid_folder+"/td.general/HHG-Spectrum.csv", dtype='str', skiprows=0)
        except:
#            print("** Error: Loading of HHG-Spectrum.csv failed.") 
            print("** Error: Loading of "+runid_folder+"/td.general/HHG-Spectrum.csv failed.")
            continue #just jump this run. Keep plotting the rest. 
        # Size of the spectrum? 
        HHG=np.asfarray(HHG)
        print("Size of HHG file: "+str(len(HHG)))

        if(len(HHG[:,0])!=HHG_samples): #NOTE: 600 is just the most employed number of samples for FFT in our computations. The size of HHG was like 600 all the time. Let's interpolate to that size if different.
            print("** Info: interpolating...")
#            HHG_func=InterpolatedUnivariateSpline(HHG[:,0], HHG[:,1])
            HHG_func=interp1d(HHG[:,0], HHG[:,1], kind='nearest') #kind='linear', 'nearest', 'zero' #WARNING: this part seems to destroy integrity of the plot

            HHG_omega_new=np.linspace(np.min(HHG[:,0]), np.max(HHG[:,0]), HHG_samples)

            Z_interp=HHG_func(HHG_omega_new)
            X=np.append(HHG_omega_new, X) #HHG frequencies for 1 run
            Z=np.append(Z_interp,      Z) #HHG amplitudes  for 1 run


            #Y=np.append(np.multiply(h*c/wavelengths1[JobID]/e,np.ones(len(HHG_omega_new))), [Y]) #we are sure that JobID is increasing monotonically. 
#            Y=np.append(np.multiply(CEP1[JobID]*180./np.pi,np.ones(len(HHG_omega_new[:,0]))),[Y])
            print("** Debug: "+str(CEP1[JobID]*180./np.pi))
            print("** Debug: checking length of new mesh: "+str(np.shape(HHG_omega_new)))
            Y=np.append(np.multiply(CEP1[JobID]*180./np.pi,np.ones(len(HHG_omega_new))),[Y])
            # HHG_samples_effective=len(HHG[:,0])
        else: #dont interpolate
            X=np.append(HHG[:,0],X)
            Z=np.append(HHG[:,1],Z)
            #Y=np.append(np.multiply(h*c/wavelengths1[JobID]/e,np.ones(len(HHG[:,0]))),[Y]) #TODO: we can use something else to plot, such as CEP or wavelengths2 !
            Y=np.append(np.multiply(CEP1[JobID]*180./np.pi,np.ones(len(HHG[:,0]))),[Y])
        print(str(int(180./np.pi * CEP1[JobID]))+" deg.")
        #print str(h*c/wavelengths1[JobID]/e)+" eV"
        del HHG

    #### ======== COMMON PART

    NumberOfRuns=len(runid)
    # Converts Y to array of reals
    Y = np.asfarray(Y)

    # Returns log10 of Z, but clean.
    logger.debug("Shape of Z: "+str(np.shape(Z)))
    Z_clipped=np.clip(Z, 0, np.max(np.max(Z))) #we remove zones lower than 0

    Z_log=np.log10(Z_clipped)
    # Z_contourf=np.reshape(Z_log, len(counter, len(Z_log)/counter))

    logger.debug("Shape of X: "+str(np.shape(X)))
    logger.debug("Shape of Y: "+str(np.shape(Y)))
    logger.debug("Number of runs: "+str(NumberOfRuns))

    Z_reshaped=np.reshape(Z_clipped, (NumberOfRuns, len(Z_clipped)//NumberOfRuns)) #integer division
    logger.debug("Shape of Z_reshaped: "+str(np.shape(Z_reshaped)))
    Z_av_CEP = np.mean(Z_reshaped, axis=0) #size matches with harmonics spectra. Then we well integrated on CEPs.
    # Z_av = np.mean(Z_clipped) #ERROR: this returns a scalar. We need a vector that depends on photon frequency.
    logger.debug("Shape of Z_clipped: " + str(np.shape(Z_clipped)))
    logger.debug("Shape of Z_av_CEP: "+str(np.shape(Z_av_CEP)))
    logger.info("CEP-average value of the HHG spectra: "+str(Z_av_CEP) )
    # We need to re-extend the list to size of Z_clipped.
    Z_av_CEP_r = np.repeat(Z_av_CEP, NumberOfRuns)
    logger.debug("Shape of repeated Z_av_CEP: "+str(np.shape(Z_av_CEP_r)))
    Z_av = np.reshape( Z_av_CEP_r, np.shape(Z_clipped))
    logger.debug("Shape of repeated Z_av_CEP: " + str(np.shape(Z_av)))
    Z_fluctuation = Z_clipped - Z_av #recommended by M. Kozak
    # DEBUG PURPOSE
    # Z = Z_av

    # THE GOOD ONE IS:
    Z = Z_fluctuation
#    Z=Z_log #Z_clipped #plotting fluctuations
    # Z=Z_clipped
    # From here, we have all the files we need. 

    xlabel=r'$\hbar \omega_{HHG}$ (eV)'
    #ylabel=r"Wavelength $\lambda$ (nm)"
    #ylabel=r"Photon energy $\hbar \omega$ (eV)"
    ylabel=r'Carrier Env. Phase (deg)'
    title=r'$HHG(\omega_{probe}$, CEP$_{pump})$'
#    cmap='Blues'
#    cmap='RdBu_r'
    # afmhot_r, hot, hot_r, bone_r, 
#    colormap='bone_r'
    colormap="RdBu_r"

    ShowFinalPlot=True
    PlotTriangles=True
    PlotGradient =False
    
    xlimit     = [2.,h*c/100e-9/e] #eV. #NOTE: we can up to 100 eV. But experiments in Si see max until 5 eV at the moment (May 2021).
    ylimit     = [-90.,90.] #deg
    numberoflevels=30
    grad_precisionX=10; grad_precisionY=10; #number of cells for interpolation for gradient computation
    Egap=2.56*e
    
    filename="HHGspectrumWithCEP-wavelength"+str(int(wavelength*1E9))+"nm-Efield"+str(int(Efield))
    
    ## EXPORT to external file.
    print("** Info: Exporting to file "+filename+".csv")
    libDatabase.ExportToTxt(np.transpose([X, Y, Z]), filename+".csv")

    #plot2dHeatMap(BrutX, BrutY, BrutZ, QuantityTitle, filename, ShowFinalPlot=True, PlotGradient=False, PlotTriangles=True, xlabel='', ylabel='', cmap_ext='afmhot_r')
    print("** Preparing plot...")
    #Signature by default: ShowFinalPlot=True, PlotGradient=False, PlotTriangles=True, xlabel='', ylabel='', cmap_ext='RdBu_r', number_of_levels=10, xlimit=[0,10], ylimit=[0,10], Egap=2.56*e, resolutionX=100, resolutionY=100, wavelength1=2656e-9
    libPlotting.plot2dHeatMap_HHG(X, Y, Z, title, filename, ShowFinalPlot, PlotGradient, PlotTriangles, xlabel, ylabel, colormap, numberoflevels, xlimit, ylimit, Egap, grad_precisionX, grad_precisionY)
    return 0
#}}}

## Selects usable data to plot HHG spectrum as function of orientation
# This routine uses pandas.
def PlotHHGsinglePulseWithOrientation(wavelength=2.65625e-6, Efield=3E9, tau=1.394500e-14, MoveIon=0,
                                      PulseType=libPlotting.PulseTypeID['SinSquared'], RiseCycles=2,
                                  InducedFields=0, LogScale=False):  # {{{
        # global Database
        #1. Find the most recent database CSV file
        # 2. Read it using pandas
        data=pd.read_csv("SiComplete-2022-01-16-BiatriSrv-v2.2.1-WithToleranceOnField-WithAll.csv", "\t", header=2)
        data.describe(include="all")
        data = data.rename(columns={"##JobID": "JobID"}) #solving hashtag problems
        data = data=data.convert_dtypes() #automatically adapts the data types for each field
        # 3. Select corresponding data
        data = data[data.Functional == libPlotting.IndexFunctional["LDA"]]  # 0: LDA, #3: TB09, #2: GGA, #1: MGGA
        data = data[data.kgrid == 40]
        data = data[data.InducedFields == 0]
        data = data.sort_values(["lambda1", "tau1", "Efield1", "Efield2", "PolarizationAngle", "kgrid", "realmesh"],
                                ascending=True)
        datatmp = data
        datatmp = datatmp[np.fabs(datatmp.lambda1 - wavelength) < 1e-9]  # Kozak 2000 nm
        datatmp = datatmp[datatmp.Efield2 == 0]
        datatmp = datatmp[datatmp.Efield1 == Efield]
        datatmp = datatmp[datatmp.CEP1 == 0]
        datatmp = datatmp[datatmp.MoveIon == MoveIon]
        datatmp = datatmp[np.fabs(datatmp.tau1 - tau)<1E-15]
        datatmp = datatmp[datatmp.STH == PulseType]
        datatmp = datatmp[datatmp.FrozenHxc == 0]
        logger.info("Number of runs available for HHG: "+str(datatmp))

        # 4. From that point, fetch ID of the run
        PolarizationAngle = datatmp.PolarizationAngle
        # logger.debug("len(RunID): "+str(len(runid)))

        # 5. Fetching the data in each subfolder's run
        X = []; Y = []; Z = []
        for JobID_id in datatmp.JobID: #{{{
            logger.debug("JobID: "+str(JobID_id))
            # logger.debug("runid: " + str(runid))
            FileLookUp = "octopus-td." + str(JobID_id) + "*.log"
            ListFound = glob.glob("*/*/*/" + FileLookUp)
            logger.debug(FileLookUp)
            runid_folder = os.path.split(os.path.abspath(ListFound[0]))[0]  # assumes only 1 finding per loop
            logger.debug(ListFound)
            logger.debug(runid_folder)

            try:
                HHG = np.loadtxt(runid_folder + "/td.general/HHG-Spectrum.csv", dtype='str', skiprows=0)
            except:
                #            print("** Error: Loading of HHG-Spectrum.csv failed.")
                logger.error("Loading of " + runid_folder + "/td.general/HHG-Spectrum.csv failed.")
                continue  # just jump this run. Keep plotting the rest.
            # Size of the spectrum?
            HHG = np.asfarray(HHG)
            logger.debug("Size of HHG file: " + str(len(HHG)))
            logger.debug("HHG_samples: "+str(HHG_samples))
            if (len(HHG[:,0]) != HHG_samples):  # NOTE: 600 is just the most employed number of samples for FFT in our computations. The size of HHG was like 600 all the time. Let's interpolate to that size if different.
                logger.info("** Info: interpolating...")
                #            HHG_func=InterpolatedUnivariateSpline(HHG[:,0], HHG[:,1])
                HHG_func = interp1d(HHG[:, 0], HHG[:, 1],
                                    kind='nearest')  # kind='linear', 'nearest', 'zero' #WARNING: this part seems to destroy integrity of the plot

                HHG_omega_new = np.linspace(np.min(HHG[:, 0]), np.max(HHG[:, 0]), HHG_samples)

                Z_interp = HHG_func(HHG_omega_new)
                X = np.append(HHG_omega_new, X)  # HHG frequencies for 1 run
                if(LogScale):
                        Z = np.append(np.log10(Z_interp+1.), Z)  # HHG amplitudes  for 1 run
                else:
                        Z = np.append(Z_interp, Z)
                # Y=np.append(np.multiply(h*c/wavelengths1[JobID]/e,np.ones(len(HHG_omega_new))), [Y]) #we are sure that JobID is increasing monotonically.
                # Y=np.append(np.multiply(CEP1[JobID]*180./np.pi,np.ones(len(HHG_omega_new[:,0]))),[Y])

                localentry = datatmp[datatmp.JobID == JobID_id]
                PolarizationAngle = np.array(localentry.PolarizationAngle)[0]
                logger.debug(PolarizationAngle * 180. / np.pi)
                logger.debug("Checking length of new mesh: " + str(np.shape(HHG_omega_new)))
                Y = np.append(np.multiply(PolarizationAngle * 180. / np.pi, np.ones(len(HHG_omega_new))), [Y])
            else:  # dont interpolate
                X = np.append(HHG[:, 0], X)
                Z = np.append(HHG[:, 1], Z)
                # Y=np.append(np.multiply(h*c/wavelengths1[JobID]/e,np.ones(len(HHG[:,0]))),[Y]) #TODO: we can use something else to plot, such as CEP or wavelengths2 !
                Y = np.append(np.multiply(PolarizationAngle * 180. / np.pi, np.ones(len(HHG[:, 0]))), [Y])
            logger.info(str(int(180. / np.pi * PolarizationAngle)) + " deg.")
            # print str(h*c/wavelengths1[JobID_id]/e)+" eV"
            del HHG
        #}}}

        #### ======== COMMON PART
        # Converts Y to array of reals
        Y = np.asfarray(Y)

        # Returns log10 of Z, but clean.
        # Z_clipped = Z #np.clip(Z, 0, np.max(np.max(Z)))  # we remove zones lower than 2
        #    Z_log=np.log10(Z_clipped)
        # Z_contourf=np.reshape(Z_log, len(counter, len(Z_log)/counter))
        # Z_av = np.mean(Z_clipped)
        # print("Average value of the distribution: " + str(Z_av))
        # Z_fluctuation = Z_clipped - Z_av  # recommended by M. Kozak
        # Z = Z_fluctuation
        #    Z=Z_log #Z_clipped #plotting fluctuations
        # Z=Z_clipped
        # From here, we have all the files we need.

        xlabel = r'$\hbar \omega_{HHG}$ (eV)'
        # ylabel=r"Wavelength $\lambda$ (nm)"
        # ylabel=r"Photon energy $\hbar \omega$ (eV)"
        ylabel = r'Sample orientation (deg)'
        title = r'$HHG(\omega_{probe}$, orientation$_{pump})$'
        #    cmap='Blues'
        #    cmap='RdBu_r'
        # afmhot_r, hot, hot_r, bone_r,
        #    colormap='bone_r'
        colormap = "RdBu_r"

        ShowFinalPlot = True
        PlotTriangles = True
        PlotGradient = False

        xlimit = [0., h * c / 200e-9 / e]  # eV. #NOTE: we can up to 100 eV. But experiments in Si see max until 5 eV at the moment (May 2021).
        ylimit = [0., 45.]  # deg
        numberoflevels = 30
        grad_precisionX = 10;
        grad_precisionY = 10;  # number of cells for interpolation for gradient computation
        Egap = 2.56 * e

        filename = "HHGspectrumWithOrientation-wavelength" + str(int(wavelength * 1E9)) + "nm-Efield" + str(int(Efield))
        if(LogScale): 
                filename += "LogScale"

        ## EXPORT to external file.
        logger.info("** Info: Exporting to file " + filename + ".csv")
        libDatabase.ExportToTxt(np.transpose([X, Y, Z]), filename + ".csv")

        # plot2dHeatMap(BrutX, BrutY, BrutZ, QuantityTitle, filename, ShowFinalPlot=True, PlotGradient=False, PlotTriangles=True, xlabel='', ylabel='', cmap_ext='afmhot_r')
        logger.info("** Preparing plot...")
        libPlotting.plot2dHeatMap_HHG(X, Y, Z, title, filename, ShowFinalPlot, PlotGradient, PlotTriangles, xlabel,
                                      ylabel, colormap, numberoflevels, xlimit, ylimit, Egap, grad_precisionX,
                                      grad_precisionY)
        return 0
#}}}

def plotHHG_batch_Efield(LogScale=False):
# HHG(omegaHHG, omega1; Efield1)
	tau=20E-15; CEP=0; MoveIon=0; InducedFields=0 #TODO: option for selecting induced fields is not included in function.
	for Efield in [5E8, 7.5E8, 1E9, 1.5E9, 2E9, 3E9, 4E9, 5E9]:
	    PlotHHGsinglePulsewithWavelength(Efield, tau, CEP, LogScale=LogScale)
    
# HHG(omegaHHG, CEP1; Efield1)
# Warning: convergence of runs should be something like 40x40x40 to be accurate.
def plotHHG_batch_CEP(LogScale=False): #{{{
    tau=20E-15; MoveIon=0; InducedFields=0; PulseType=3; RiseCycles=2
    for wavelength in [800e-9, 1600e-9, 2200E-9, 3200E-9]:
        for Efield in [5E8, 1E9, 2E9, 3E9, 7.5E9]:
            try: 
                PlotHHGsinglePulsewithCEP(wavelength, Efield, tau, MoveIon, PulseType, RiseCycles, InducedFields)
            except:
                logger.warning("No data for (wavelength, field)="+str(wavelength*1E0)+"nm, "+str(Efield*1E-9)+"V/nm.")
                continue
#}}}

# plotHHG_batch_CEP() #NOTE: not sure there is an effect. 200'000 corehours submitted to Salomon on Oct. 22nd 2019.

## HHG spectra computations as function of CEP for Kozak et al 2022.
# plotHHG_batch_Efield() #Initially related to Derrien et al., https://arxiv.org/abs/2104.08971
def KozakOptica2022_batch_CEP(): #{{{
    for Efield in [3E9, 1.05E10]:
        PlotHHGsinglePulsewithCEP(wavelength=2656e-9, Efield=Efield, tau=13.945E-15, MoveIon=0,
                                  PulseType=libPlotting.PulseTypeID['SinSquared'], RiseCycles=0, InducedFields=0) #HHG spectra
    #    PlotHHGsinglePulsewithCEP(2656e-9, Efield, 13.945E-15, 0, 1, 0, 0)
#}}}

def KozakOptica2022_batch_orientation(LogScale=False): #{{{
    # for Efield in [3E9, 1.05E10]:
    for Efield in [9E9]:
        PlotHHGsinglePulseWithOrientation(wavelength=2656e-9, Efield=Efield, LogScale=LogScale)

## HHG spectra computation as function of angle of incidence (Kozak, Derrien et al 2022).
if(__name__=="__main__"):
    LogScale=True
    # KozakOptica2022_batch_orientation(LogScale=LogScale)
    KozakOptica2022_batch_CEP() #no logscale possible here (we look at fluctuations)
