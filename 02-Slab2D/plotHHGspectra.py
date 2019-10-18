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

X=[]; Y=[]; Z=[]
counter=0

for FolderName in glob.glob(path): #{{{
    print Header+"** Info: entering folder: "+str(FolderName)+"\n"
    # If the run considers only single pulses
    filename = FolderName+"../inp"
    Efield2_SI = libPlotting.GetParameterValueFromInputFile(filename, "EfieldBsi = ")
    PulseDuration_SI = libPlotting.GetParameterValueFromInputFile(filename, "tauAsi = ")

    print "Pulse duration 1 (s): "+str(PulseDuration_SI)
    print "Field amplitude 2 (V/m): "+str(Efield2_SI)

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

# From here, we have all the files we need. 
ShowFinalPlot=False
libPlotting.plot2dHeatMap(X, Y, Z, r"HHG(\omega_{\mbox{probe}}, \lambda_{\mbox{pump}})", "HHGspectrum", ShowFinalPlot, False)
