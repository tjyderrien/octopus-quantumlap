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

## @package plotHHGspectra
# Collects the generated time-dependent octopus simulations spectra HHG.dat from all subfolders (can last minutes)
# Then prepares a picture using the collected data for each identified field amplitude

import glob
import os

from spp_extended_theory.Keldysh.libKeldysh import *
from scipy.interpolate import interp1d

from Libs import libPlotting

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':'14'})
## for Palatino and other serif fonts use:
#rc('font', **{'family':'serif', 'serif':['Palatino'], 'size':'16'})
rc('text', usetex=False) #True: Does not work on Draco. Just put False, then, and it will be smooth. 
mp.rcParams['legend.numpoints'] = 1

Header="[plotHHGspectra] "

ShowFinalPlot              = 1
HHG_samples                = 600 #*2, *3, *4

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
	        HHG=loadtxt(FolderName+"../td.general/HHG-Spectrum.csv", dtype='str', skiprows=1)
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
    

Database = np.loadtxt(DatabaseFile, skiprows=0)
## ======== METHOD 2: construct the sub-database of interest

def PlotHHGsinglePulsewithWavelength(Efield, tau=20E-15, CEP=0, MoveIon=0, PulseType=3, RiseCycles=2): #{{{
    
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
            HHG=loadtxt(runid_folder+"/td.general/HHG-Spectrum.csv", dtype='str', skiprows=1)
        except:
            print("** Error: Loading of HHG-Spectrum.csv failed.") 
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
    Z_log=np.log10(Z_clipped)
    # Z_contourf=np.reshape(Z_log, len(counter, len(Z_log)/counter))
    Z=Z_log #Z_clipped
    # Z=Z_clipped

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

    ## EXPORT to external file.
    print("** Info: Exporting to file "+filename+".csv")
    ExportToTxt(np.transpose([X, Y, Z]), filename+".csv")

    print("** Preparing plot...")
    libPlotting.plot2dHeatMap_HHG(X, Y, Z, title, filename, ShowFinalPlot, PlotGradient, PlotTriangles, xlabel, ylabel, colormap, numberoflevels, xlimit, ylimit, Egap)
    return 0
#}}}

## Prepares color plot HHG(omega_HHG, CEP) 
def PlotHHGsinglePulsewithCEP(wavelength=800e-9, Efield=1E9, tau=20E-15, MoveIon=0, PulseType=3, RiseCycles=2, InducedFields=0): #{{{
    
    # Declaring variables
    X=[]; Y=[]; Z=[]
   
    # FilterSinglePulseDatabaseWithParameters_wavelengths(Database, Efield1, pulseDuration, CEP, PulseType=3, RiseTime=2)
    #Efield =7.5E8 #V/m
    #tau    =20e-15
    #CEP    =0
    #MoveIon=0
    #PulseType=3; RiseCycles=2
    Filtered1= libPlotting.FilterSinglePulseDatabaseWithParameters_CEP(Database, wavelength, Efield, tau, PulseType, RiseCycles, InducedFields)
    Filtered1=np.asfarray(Filtered1) #lets use only numbers there.
    # SORTING
    # One must sort here the DB as function of the sequence we would like to plot.
    # Sorted1=np.sort(Filtered1, axis=0)
    Sorted2=Filtered1[np.lexsort((Filtered1[:, libPlotting.IndexTable['Efield1']], Filtered1[:,
                                                                                   libPlotting.IndexTable['CEP1']]))]
    # NOTE: how to manage the duplicated computations? We would like to keep longest computations. But may not remove duplicates. For now: manually disabling duplicated entries in DB. 
    # One could use np.unique() function, but we may loose the JobID then? 
    # Sorted3=libDatabase.FilterDatabaseRemove( #lets remove the shortest duplicated run

    # Extracting the data from filtered database
    try: 
        runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, Energy_max, Wpi, PulseType, RiseTime, FinalTime  = libPlotting.ExtractData(Sorted2)
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
        FileLookUp="octopus-td."+str(runid_int)+"*.log"
        ListFound=glob.glob("*/*/"+FileLookUp)
        print(FileLookUp)
        runid_folder=os.path.split(os.path.abspath(ListFound[0]))[0] #assumes only 1 finding per loop
        print(runid_folder)
        try:
            HHG=loadtxt(runid_folder+"/td.general/HHG-Spectrum.csv", dtype='str', skiprows=1)
        except:
            print("** Error: Loading of HHG-Spectrum.csv failed.") 
            continue #just jump this run. Keep plotting the rest. 
        # Size of the spectrum? 
        HHG=np.asfarray(HHG)
        print("Size of HHG file: "+str(len(HHG)))

        if(len(HHG[:,0])!=HHG_samples): #NOTE: 600 is just the most employed number of samples for FFT in our computations. The size of HHG was like 600 all the time. Let's interpolate to that size if different.
            print("** Info: interpolating...")
#            HHG_func=InterpolatedUnivariateSpline(HHG[:,0], HHG[:,1])
            HHG_func=interp1d(HHG[:,0], HHG[:,1], kind='nearest') #kind='linear', 'nearest', 'zero'

            HHG_omega_new=np.linspace(np.min(HHG[:,0]), np.max(HHG[:,0]), HHG_samples)

            Z_interp=HHG_func(HHG_omega_new)
            X=np.append(HHG_omega_new, X) #HHG frequencies for 1 run
            Z=np.append(Z_interp,      Z) #HHG amplitudes  for 1 run


            #Y=np.append(np.multiply(h*c/wavelengths1[JobID]/e,np.ones(len(HHG_omega_new))), [Y]) #we are sure that JobID is increasing monotonically. 
#            Y=np.append(np.multiply(CEP1[JobID]*180./np.pi,np.ones(len(HHG_omega_new[:,0]))),[Y])
            print("** Debug: "+str(CEP1[JobID]*180./np.pi))
            print("** Debug: checking length of new mesh: "+str(np.shape(HHG_omega_new)))
            Y=np.append(np.multiply(CEP1[JobID]*180./np.pi,np.ones(len(HHG_omega_new))),[Y])
        else: #dont interpolate
            X=np.append(HHG[:,0],X)
            Z=np.append(HHG[:,1],Z)
            #Y=np.append(np.multiply(h*c/wavelengths1[JobID]/e,np.ones(len(HHG[:,0]))),[Y]) #TODO: we can use something else to plot, such as CEP or wavelengths2 !
            Y=np.append(np.multiply(CEP1[JobID]*180./np.pi,np.ones(len(HHG[:,0]))),[Y])
        print(str(int(180./np.pi * CEP1[JobID]))+" deg.")
        #print str(h*c/wavelengths1[JobID]/e)+" eV"
        del HHG

    #### ======== COMMON PART 
    # Converts Y to array of reals
    Y = np.asfarray(Y)

    # Returns log10 of Z, but clean.
    Z_clipped=np.clip(Z, 2, np.max(Z))
    Z_log=np.log10(Z_clipped)
    # Z_contourf=np.reshape(Z_log, len(counter, len(Z_log)/counter))
    Z=Z_log #Z_clipped
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
    colormap='bone_r'

    ShowFinalPlot=False
    PlotTriangles=False
    PlotGradient =False
    
    xlimit     = [0,50] #eV
    ylimit     = [-90,90] #deg
    numberoflevels=30
    Egap=2.56*e
    
    filename="HHGspectrumWithCEP-wavelength"+str(int(wavelength*1E9))+"nm-Efield"+str(int(Efield))
    
    ## EXPORT to external file.
    print("** Info: Exporting to file "+filename+".csv")
    ExportToTxt(np.transpose([X, Y, Z]), filename+".csv")

    #plot2dHeatMap(BrutX, BrutY, BrutZ, QuantityTitle, filename, ShowFinalPlot=True, PlotGradient=False, PlotTriangles=True, xlabel='', ylabel='', cmap_ext='afmhot_r')
    print("** Preparing plot...")
    libPlotting.plot2dHeatMap_HHG(X, Y, Z, title, filename, ShowFinalPlot, PlotGradient, PlotTriangles, xlabel, ylabel, colormap, numberoflevels, xlimit, ylimit, Egap)
    return 0
#}}}

def plotHHG_batch_Efield():
# HHG(omegaHHG, omega1; Efield1)
	tau=20E-15; CEP=0; MoveIon=0; InducedFields=0 #TODO: option for selecting induced fields is not included in function.
	for Efield in [5E8, 7.5E8, 1E9, 1.5E9, 2E9, 3E9, 4E9, 5E9]:
	    PlotHHGsinglePulsewithWavelength(Efield, tau, CEP)
    
# HHG(omegaHHG, CEP1; Efield1)
def plotHHG_batch_CEP():
	tau=20E-15; MoveIon=0; InducedFields=0; PulseType=3; RiseCycles=2
        for wavelength in [800e-9, 1600e-9, 2200E-9, 3200E-9]:
            for Efield in [5E8, 1E9, 2E9, 3E9, 7.5E9]:
                try: 
                    PlotHHGsinglePulsewithCEP(wavelength, Efield, tau, MoveIon, PulseType, RiseCycles, InducedFields)
                except:
                    print(Header+"** Warning: no data for (wavelength, field)="+str(wavelength*1E0)+"nm, "+str(Efield*1E-9)+"V/nm.")
                    continue

# plotHHG_batch_CEP() #NOTE: not sure there is an effect. 200'000 corehours submitted to Salomon on Oct. 22nd 2019.

plotHHG_batch_Efield() 

