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

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':'14'})
## for Palatino and other serif fonts use:
#rc('font', **{'family':'serif', 'serif':['Palatino'], 'size':'16'})
rc('text', usetex=True) #True: Does not work on Draco. Just put False, then, and it will be smooth. 
mp.rcParams['legend.numpoints'] = 1

Header="[plotHHGspectra] "

ShowFinalPlot              = 0
path = "./"
for FolderName in glob.glob(path):
    print Header+"** Info: entering folder: "+str(FolderName)+"\n"
    # If the run considers only single pulses
    filename = path+"../inp"
    Efield2_SI = libPlotting.GetParameterValueFromInputFile(filename, "EfieldBsi = ")
    
    if(Efield2_SI < 1): #we only select runs that are with single pulses. 
        # NOTE: We need multiple wavelengths, but one spectrum per wavelength. 
    
        # List all the files HHG-spectrums.csv in subfolders
        HHG=loadtxt(path+"../td.general/HHG-spectrum.csv", dtyoe='str', skiprows=1)
        
        # Capture the spectrum
        print HHG
        
        # Interpolate it to the size of the others
        # Add it as a line in table
        # Plot the final picture
        # Add analytical Stark effect there
