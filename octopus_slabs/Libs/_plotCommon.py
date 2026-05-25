#!/usr/bin/env python
#-*- coding: utf-8 -*-

# Copyright (C) 2016-2023 T. J.-Y. Derrien
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

## @module _plotCommon
# Shared imports, constants, matplotlib setup, and database initialisation
# used by all plotFinalQuantities sub-modules.

## IMPORT LIBRARIES
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline
#import matplotlib.pyplot as plt
# from pylab import *

import matplotlib  as mp
from   matplotlib import rc
from   mpl_toolkits.axes_grid1.inset_locator import inset_axes
# from matplotlib.tri import Triangulation, UniformTriRefiner
from scipy.optimize import curve_fit
from scipy.constants import h, c, e, hbar, m_e, epsilon_0, pi
import numpy as np
import matplotlib.pyplot as plt

import os
# import pdb
from octopus_slabs.Libs import libExcitationDatabase, libString, libPlotting
# import libKeldysh as Keldysh
# from libAtomicUnits import *
# from libPlottingSettings import *
from octopus_slabs.Libs import libPlotting
from octopus_slabs.Libs.libPlotting import IndexTable, IndexFunctional
from octopus_slabs.Libs.libLogging import init_logger

# libKeldysh = import_module("spp_extended_theory.Keldysh.libKeldysh")
from spp_extended_theory.Keldysh import libKeldysh, libKeldyshZhukov
# from libKeldysh import *
from spp_extended_theory.Libs import libDatabase
from spp_extended_theory.Libs import libImportOpticalData as libOpticalData

#logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
#logging.basicConfig(level=logging.INFO, handlers=[logging.StreamHandler()]) #filename='plotFinalQuantities.log', 
logger = init_logger(__name__, verbose=False) #"plotFinalQuantities")


rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':'14'})
## for Palatino and other serif fonts use:
#rc('font', **{'family':'serif', 'serif':['Palatino'], 'size':'16'})
rc('text', usetex=True) #True: Does not work on Draco. Just put False, then, and it will be smooth. 
mp.rcParams['legend.numpoints'] = 1

Header="[plotFinalQuantities] "

# EnablePeakFinder           = 0 #TODO: base it on presence of file in code root. Easier for configuration.
myhost = os.uname()[1]
logger.info(myhost)
if("biatri-srv" in myhost):
    EnablePeakFinder       = 0 #valid for BIATRI-SRV
else: 
    EnablePeakFinder       = 1
ShowFinalPlot              = 0

BicolorE1E2plotDictionnary={'Nexc': 0, 'E2': 4, 'J': 1, 'Energy': 2, 'Wpi': 3, 'NexcGrad': 5}

Density=5E28 #specific to Si #TODO: compute automatically? Given by pseudo-potential...
n0 = 4.*Density #maximum density of electrons in Si
xscale = 1E-4 #convert W/m2 to W/cm2.
try:
    SampleVolume = libString.GetCellVolume()
except:
    logger.warning((Header + "** Warning: sample computation volume was failed!"))
    SampleVolume = 0

## These data should be moved
Egap_SI = 2.56*e #LDA: 2.56 eV, TB09: 3.04 eV. [Waroquiers and Gonze, Phys. Rev. B 87, 075121 (2013)] #TODO: this will depend on which functional we have used! TB09 or LDA !
meff_Waroquiers_LDA = 0.2226 #0.2226 #Waroquiers, D. et al. Physical Review B, 2013, 87, 075121.
meff = meff_Waroquiers_LDA

# Capture the data
# Several locations are possible.
path3='*/*/*/td.general/'
path2='*/*/td.general/'
path1='*/td.general/'
path0="./td.general/"

UseThisPath=path3

if(EnablePeakFinder): 
  from scipy.signal import find_peaks

PlottingTable = libExcitationDatabase.GetExcitationDatabase(UseThisPath) #completes the existing database, or fetches the new one in all the subfolders where this is executed. 

logger.info(PlottingTable)

#PlottingTable = libDatabase.FilterDatabase(PlottingTable, 0, IndexTable['FieldInduced']) #removes all induced fields data

#: TODO one has to extract the volume! But how to proceed when it is time-dependent? 

runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingTable)

logger.info(Header+"** Preparing the summary plots...")

#### =============== PLOTTING STARTS HERE ==================
# The filtering will be done on the fly. Then we are sure of 
# what we plot in each figure. 

### Defining plotting structure

# Observables:        {N_exc^max, N_exc^final, irreversibility, TotalCurrent_final, Energy_max, Energy_final}

#### Plotting: parameters for controlling the canvas
# if(EnableDoublePulseAnalysis==0):
  # 2 plots
  
output_format='publication' #'slides' #'publication'

plotA=211; plotB=212; #plotC=313
# elif (wavelength_numbers > 1):
# plotA=311; plotB=312; plotC=313;
# else: #if(EnableDoublePulseAnalysis==1):
# 3 plots
# plotA=311; plotB=312; plotC=313;
# 4 plots 
# plotA=221; plotB=222; plotC=223; plotD=224;

# Size of main canvas
GoldNumber = ((1.+np.sqrt(5.))/2.)
SizeX = 12.
#SizeY = SizeX / GoldNumber
SizeY = 6.
OneSizeFitsAll="small"  #caption size for all legends of physics or rocknroll

PulseDurationLocal = 0e0

