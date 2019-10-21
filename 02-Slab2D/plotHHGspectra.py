#!/usr/bin/env python2
#-*- coding: utf-8 -*-

# Copyright (C) 2016-2019 T. J.-Y. Derrien
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

# @module plotHHGspectra
# Collects the generated time-dependent octopus simulations spectra HHG.dat from all subfolders (can last minutes)
# Then prepares a picture using the collected data for each identified field amplitude

from libKeldysh import *
import libPlotting
import glob

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':'14'})
## for Palatino and other serif fonts use:
#rc('font', **{'family':'serif', 'serif':['Palatino'], 'size':'16'})
rc('text', usetex=False) #True: Does not work on Draco. Just put False, then, and it will be smooth. 
mp.rcParams['legend.numpoints'] = 1

Header="[plotHHGspectra] "

ShowFinalPlot              = 1

path = "*/td.general/"
DatabaseFile="Database.csv"

X=[]; Y=[]; Z=[]
counter=0

# NOTE: two ways to capture the corresponding data can be employed. 
#
# 1. Looking the individual input file.
# 2. Or accessing the database using the ID of the run. Must be much more clever way.
    
## NOTE: METHOD 1: accessing each subfolder. 
def ParseAllFolders():
	for FolderName in glob.glob(path): #{{{ #NOTE: this takes everything available: not sorted.
	    print Header+"** Info: entering folder: "+str(FolderName)+"\n"
	  
	    ### METHOD 1: rebuilding database without choosing the order of the sequence.  
	    filename = FolderName+"../inp"
	    Efield1_SI = libPlotting.GetParameterValueFromInputFile(filename, "EfieldAsi = ")
	    Efield2_SI = libPlotting.GetParameterValueFromInputFile(filename, "EfieldBsi = ")
	    PulseDuration_SI = libPlotting.GetParameterValueFromInputFile(filename, "tauAsi = ")
	    PulseDuration2_SI= libPlotting.GetParameterValueFromInputFile(filename, "tauBsi = ")
	    wavelength1_SI   = libPlotting.GetParameterValueFromInputFile(filename, "lambdaA = ")
	    wavelength2_SI   = libPlotting.GetParameterValueFromInputFile(filename, "lambdaB = ") 
   
	    print "Pulse duration 1 (s): "+str(PulseDuration_SI)
	    print "Field amplitude 2 (V/m): "+str(Efield2_SI)
	
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
	        print "Size of HHG file: "+str(len(HHG))
	#        print HHG
	#        print "Extracting 1 column: "
	#        print HHG[:,1]
	        # Interpolate it to the size of the others
	        
	        # Add it as a line in table
	        X=np.append(HHG[:,0], X) #HHG frequencies for 1 run
	        Z=np.append(HHG[:,1], Z) #HHG amplitudes  for 1 run
	        print counter
	        Y=np.append(counter*np.ones(len(HHG)), [Y]) #Y should contain the parameter we variate. If taking it all, we use a counter to represent them all as a simple stack for now.
	        print Y
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
	        print "** Info: Double pulses are ignored."
	# }}} #glob loop

#### ======== METHOD 2: construct the sub-database of interest

Database = np.loadtxt(DatabaseFile, skiprows=0)
# FilterSinglePulseDatabaseWithParameters(Database, Efield1, pulseDuration, CEP, PulseType=3, RiseTime=2)
Efield1=1E9 #V/m
tau1   =30e-15
CEP1   =0
Filtered1=libPlotting.FilterSinglePulseDatabaseWithParameters(Database, Efield1, tau1, CEP1)

# One can sort here the DB as function of the sequence we would like to plot.
# Filtered1 

# Extracting the data from filtered database
runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, Wpi, PulseType, RiseTime, FinalTime  = libPlotting.ExtractData(Filtered1)

for JobID in np.arange(0,len(runid)):
    # Lookup for the corresponding file. Reconstruct the folder name: faster than finding the file in TB or data? 
    RebuiltFolderName="lambda"+str(wavelengths1[JobID])
    RebuiltFolderName+="-tau"+str(pulses1[JobID])
    RebuiltFolderName+="-field"+str(Efield1[JobID])
    RebuiltFolderName+="-CEP"+str(CEP1[JobID])
    RebuiltFolderName+="-"+str(PulseTypeID_rev[PulseType[JobID]])
    RebuiltFolderName+=str(RiseTime)
    RebuiltFolderName+=str(FreezeHXC)+"/"
    print RebuiltFolderName
#    RebuiltFolderName+=""lambda"${lambda2}-"tau"${tau2}-"field"${Efield2}-"CEP"${CEP2}-"delay"${delay}-${PulseShape}${CrossPolarizationText}${GaugeFieldText}${MoveIonsChain} #NOTE appears a bit complicated. 

    # HHG=loadtxt(FolderName+"../td.general/HHG-Spectrum.csv", dtype='str', skiprows=1)
                # except:
                #  print "** Error: Loading of HHG-Spectrum.csv failed." 
                # Size of the spectrum? 
    #            HHG=np.asfarray(HHG)
        #        print "Size of HHG file: "+str(len(HHG))
        #        print HHG
        #        print "Extracting 1 column: "
        #        print HHG[:,1]
                # Interpolate it to the size of the others

                # Add it as a line in table
      #          X=np.append(HHG[:,0], X) #HHG frequencies for 1 run
      #          Z=np.append(HHG[:,1], Z) #HHG amplitudes  for 1 run
      #          print counter
      #          Y=np.append(counter*np.ones(len(HHG)), [Y]) #Y should contain the parameter we variate. If taking it all, we use a counter to represent them all as a simple stack for now.
      #          print Y
                # Prepare the final table
                # Plot the intermediate picture
                # plt.figure()
                # plt.plot(X, Z)
                # plt.show() 
                # Add analytical Stark effect there
                # Cleaning
      #          del HHG
      #          counter+=1


# Once has the data, let's extract the job numbers, extract the HHG-Spectrum.csv files and build the figure.
# libDatabase.ExtractMaterialData(DataseFile)
# libDatabase.ExtractDataDb(DatabaseFile)

exit()

#### ======== COMMON PART 

# Returns log10 of Z, but clean.
Z_clipped=np.clip(Z, 1, np.max(Z))
Z_log=np.log10(Z_clipped)
Z_contourf=np.reshape(Z_log, (counter, len(Z_log)/counter))
Z=Z_clipped
# From here, we have all the files we need. 
xlabel=r"$\hbar \omega_{HHG}$ (eV)"
ylabel="Runs"
ShowFinalPlot=False
PlotTriangles=False
plt.figure()
# plt.contourf(X,Y,Z_contourf)
plt.savefig("HHG.png")
# plt.show()
libPlotting.plot2dHeatMap(X, Y, Z, r"$HHG(\omega_{probe}$, $\lambda_{pump})$", "HHGspectrum", ShowFinalPlot, False, PlotTriangles, xlabel, ylabel)
