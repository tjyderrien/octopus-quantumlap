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

## @module plotPaperCollections
# High-level convenience functions that assemble publication-ready figures
# by orchestrating calls to the lower-level plotting modules.
# Functions: plot_Zhukov_SingleColor_SiO2, plot_Zhukov_SingleColor_LDA,
#            plot_Zhukov_SingleColor_TB09, plot_Zhukov_BicolorMaps,
#            plot_ShiftCurrents, plot_HHG, plot_Vorobyev_Au

from octopus_slabs.Libs._plotCommon import *  # noqa: F401,F403
from octopus_slabs.Libs.plotSinglePulse import (
    plotEffectOfRiseTime,
    plotEffectOfIntensity_wavelengths,
    plotEffectOfIntensity_PulseDuration,
    plotEffectOfIntensity_PulseDuration_Metals,
    plotEffectOfIntensity_CEP,
    plotEffectOfCEP,
)
from octopus_slabs.Libs.plotWavelength import (
    plotEffectOfWavelength_Tolerance,
    Bicolor_plotEffectOfWavelength_Tolerance,
    Bicolor_plotEffectOfCEP,
    plotEffectOfWavelength_2dMap,
)
from octopus_slabs.Libs.plotDoublePulse import (
    plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple,
    plotAdjustedCollisionFrequencies_wavelength,
)

def plot_Zhukov_SingleColor_SiO2(): #{{{
    Egap_SI = 9.*e; Egap_SI_indirect=Egap_SI #we don't care of indirect gap
    PulseType  = 3 #STH: 2, StH: 3
    RiseCycles = 2
    ShowPlot = 1
    NumberOfAtoms=9. #3 Si & 6 Oxygen. 
    plotEffectOfWavelength_Tolerance(PlottingTable, FieldTolerance=1E-4, pulseDuration=20e-15, CEP=0.,
                                     PulseType=PulseType, RiseTimeRequested=RiseCycles, ShowFinalPlot=ShowPlot,
                                     NumberOfAtoms=NumberOfAtoms, Egap_SI=Egap_SI, Egap_SI_indirect=Egap_SI_indirect,
                                     wavelengths1_min=500e-9, wavelengths1_max=3300e-9)

    #def plotEffectOfWavelength_Tolerance(PlottingTable, FieldTolerance=0.05, pulseDuration=20e-15, CEP=0., PulseType=3,
    #                                     RiseTimeRequested=2, ShowFinalPlot=0, NumberOfAtoms=2., Egap_SI=2.56 * e,
    #                                     Egap_SI_indirect=0.51 * e, wavelengths1_min=200e-9, wavelengths1_max=3300e-9,
    #                                     ShowDensities=False, PlotIrreversibility=False, PlotFreeCarrierAbsorption=True,
    #                                     EnergyAxis=True, Functional="LDA"):  # {{{
#}}}                                            

## Prepare all plots for single-color excitation with VP Zhukov [Phys. Rev. B / PRL]
def plot_Zhukov_SingleColor_LDA(EnableMultiwavelengthPlot=False): #{{{
    Egap_SI = 2.56*e
    #Egap_SI = 1.55*e #Dressed in presence of driving photon field
    #Egap_SI = 1.5*e #Manually adjusted gap. 
    Egap_SI_indirect = 0.51*e
    PulseType  = 3 #STH: 2, StH: 3
    wavelength = 800e-9
    ShowPlot   = 0
    PlotWpiAxis = True
    MoveIon=0
    NumberOfAtoms=2. #Si
    
    ## These are experimental data from Palik. More correct would be to use linear response computed from TDDFT. 
    OpticalIndex={ #
        483e-9: 4.4028, 484e-9: 4.3964,
        800e-9: 3.6924167231905, 1030e-9: 3.61279196759238,
                 1600e-9: 3.4826,
        2200e-9: 3.4548, 2600e-9: 3.4457, 3200e-9: 3.43438589865999}
    
    ## NOTE: DUMMY OpticalIndex account
    #OpticalIndex={
        #483e-9: 1., 484e-9: 1.,
        #800e-9: 1., 1600e-9: 1., 
        #2200e-9: 1., 2600e-9: 1., 3200e-9: 1. }
    
    #Egap_SI_dressed = {483e-9:, 484e-9: ,
                       #800e-9: 1.55*e, 1600e-9: ,
                       #2200e-9:, 2600e-9: }
    
    ## Plots effect of the STH Rise Time. We even found an incomplete computations thanks to this. Works super well, then. 
    RiseCycles = 0 #0: undefined, mixes everything. Basically, we must redefine all runs which had no well defined RiseCycle so far. 
    Efield     = 1E9
    #plotEffectOfRiseTime(PlottingTable, wavelength, 20e-15, Efield, PulseType, ShowPlot, Egap_SI)
    
    # This is a manually prepared figure. 

    #plotAdjustedCollisionFrequencies_wavelength()
    
    ## Effects of CEP on single pulse excitation, current and energy, and harmonic spectrum
    RiseCycles = 2
    Efield     = 1E9
    # plotEffectOfCEP(PlottingTable, wavelength, 20e-15, Efield, PulseType, RiseCycles, ShowPlot, Egap_SI)
    #plotEffectOfCEP(PlottingTable, wavelength, 10e-15, Efield, PulseType, RiseCycles, ShowPlot)
    Efield     = 5E9 #NOTE: Saturation may prevent to see any difference. 
    # plotEffectOfCEP(PlottingTable, wavelength, 20e-15, Efield, PulseType, RiseCycles, ShowPlot, Egap_SI)
    #plotEffectOfCEP(PlottingTable, wavelength, 10e-15, Efield, PulseType, RiseCycles, ShowPlot)
    #NOTE: === LOCAL FILTERING
    #PlottingTable_filtered = FilterDatabaseRemove(PlottingTable, 30e-15, IndexTable['pulses1'])
    PlottingTable_filtered = PlottingTable
    

    ## Effect of intensity figures - not much interesting for "long" pulses
    # plotEffectOfIntensity_CEP(PlottingTable_filtered, wavelength, 20e-15, False, 0.0, PulseType, RiseCycles, ShowPlot, Egap_SI)
    #plotEffectOfIntensity_CEP(PlottingTable_filtered, wavelength, 10e-15, False, 0.0, PulseType, RiseCycles, ShowPlot)
    #plotEffectOfIntensity_CEP(PlottingTable_filtered, wavelength, 30e-15, False, 0.0, PulseType, RiseCycles, ShowPlot)
    
    ### ABSORPTION PROBABILITY FOR 800 nm - Phys Rev B Derrien et al 2020
    #plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI, Egap_SI_indirect, OpticalIndex[wavelength]) #This one construct the ionization probabilities! 
    
    ## RESONANT ABSORPTION - not published, not understood
    wavelength = 483e-9 #Resonant wavelength results for LDA.
    #plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI, Egap_SI_indirect, OpticalIndex[wavelength]) #This one construct the ionization probabilities! 
    wavelength = 484e-9 #Resonant wavelength results for LDA.
    #plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI, Egap_SI_indirect, OpticalIndex[wavelength]) #This one construct the ionization probabilities! 
    
    ## READY
    wavelength = 800e-9
    #plotEffectOfIntensity_wavelengths(PlottingTable_filtered, 20e-15, 0., False, wavelength, PulseType, RiseCycles, ShowPlot, Egap_SI) 
    ##plotEffectOfIntensity_wavelengths(PlottingTable_filtered, 10e-15, 0., False, wavelength, PulseType, RiseCycles, ShowPlot)
    #plotEffectOfIntensity_wavelengths(PlottingTable_filtered, 30e-15, 0., False, wavelength, PulseType, RiseCycles, ShowPlot, Egap_SI)
    
    ## CHECKING DRUDE MODEL VALIDITY (results as function of wavelength)
    ShowPlot=1
    #plotEffectOfWavelength_2dMap(PlottingTable, 20e-15, PulseType, RiseCycles, ShowPlot)
    
    ### NOTE THIS IS FIG3 OF THE PRL
    if(EnableMultiwavelengthPlot):
        plotEffectOfWavelength_Tolerance(PlottingTable, FieldTolerance=1E-4, pulseDuration=20e-15, CEP=0.,
                                         PulseType=PulseType, RiseTimeRequested=RiseCycles, ShowFinalPlot=ShowPlot,
                                         NumberOfAtoms=NumberOfAtoms, Egap_SI=Egap_SI,
                                         Egap_SI_indirect=Egap_SI_indirect, wavelengths1_min=500e-9,
                                         wavelengths1_max=3300e-9)
        # def plotEffectOfWavelength_Tolerance(PlottingTable, FieldTolerance=0.05, pulseDuration=20e-15, CEP=0., PulseType=3,
        #                                     RiseTimeRequested=2, ShowFinalPlot=0, NumberOfAtoms=2., Egap_SI=2.56 * e,
        #                                     Egap_SI_indirect=0.51 * e, wavelengths1_min=200e-9, wavelengths1_max=3300e-9,
        #                                     ShowDensities=False, PlotIrreversibility=False, PlotFreeCarrierAbsorption=True,
        #                                     EnergyAxis=True, Functional="LDA"):  # {{{
    
    ## Checking what is happening in the mid-IR world. NOTE: NOT ENOUGH DATA? 
    #plotEffectOfWavelength_Tolerance(PlottingTable, 1E-4, 60e-15, 0., PulseType, RiseCycles, ShowPlot, Egap_SI, Egap_SI_indirect, 3e-6, 25e-6)
    #plotEffectOfWavelength_Tolerance(PlottingTable, 0.05, 30e-15, 0., PulseType, RiseCycles, ShowPlot, Egap_SI, Egap_SI_indirect) #ALMOST NO DATA for this one

    ### OLD STUFFS
    #plotEffectOfIntensity_CEP(PlottingTable, wavelength, 10e-15, PulseType, RiseCycles, ShowPlot)
    ## Plots the N_exc(wavelength, field) (2D plot) # NOTE: this can be compared with GrojoPRB2013. 
    
    ## MID-IR: Preparation of the mid-IR experiment: measuring the cross section absorption in mid-IR 
    #wavelength = 1600e-9
    #plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, 1600e-9, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI, Egap_SI_indirect, OpticalIndex[wavelength]) #contains MoveIon, yes?
    
    # Some data are available for 2200 nm, but they are very sparse (not priority)
    #plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, 2200e-9, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI, Egap_SI_indirect, OpticalIndex[wavelength]) #contains MoveIon, yes? 
    
    #plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, 2600e-9, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI, Egap_SI_indirect, OpticalIndex[wavelength]) #contains MoveIon, yes?
    #plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, 3200e-9, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI, Egap_SI_indirect) #contains MoveIon, yes? 
    
    #for wavelength in [3200E-9, 2600e-9, 2200e-9, 1600e-9, 800e-9]: 
    # for wavelength in [3200e-9, 1600e-9, 800e-9, 483e-9]: #, 1600e-9, 3200E-9]: #, 1600e-9, 800e-9]:
    for wavelength in [1030e-9]:  # , 1600e-9, 3200E-9]: #, 1600e-9, 800e-9]:
      #try:
      plotEffectOfIntensity_PulseDuration(PlottingTable_filtered,  wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI, Egap_SI_indirect, OpticalIndex[wavelength], MoveIon, "LDA")
      #except:
        #print Header+"** Warning: no data could be found, or bug took place for "+str(wavelength*1E9)+" nm."
    
    material="Si"
    PrimitiveCellArea = {'Si': 1E-20} #WARNING: very approximative!!!
    RatioExternal     = {'Si': 1.}    #Useful only for induced fields cases
    ## Plotting absorbed energy also
    for PulseType in [libPlotting.PulseTypeID['SinSquared'], libPlotting.PulseTypeID['StH']]:
        # for wavelength in [3200e-9, 1600e-9, 1030e-9, 800e-9]:
        for wavelength in [1030E-9]:
        #0: Gaussian; 1: SinSquared
            try: 
                plotEffectOfIntensity_PulseDuration_Metals(PlottingTable_filtered,  wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, PrimitiveCellArea[material], OpticalIndex[wavelength], RatioExternal[material], MoveIon)
            except: 
                logger.warning((Header+"Plotting of Energy for Si failed for wavelength "+str(int(1E9*wavelength))+", pulsetype: "+str(PulseType)))
    
#}}}

## Plotting excitation rate as function of intensity    
def plot_Zhukov_SingleColor_TB09(): #{{{                                            
    OpticalIndex={ #
        407.58e-9: 5.3705, 483e-9: 4.4028, 484e-9: 4.3964,
        800e-9: 3.6924167231905, 1030e-9: 3.61279196759238, 1600e-9: 3.4826,
        2200e-9: 3.4548, 2600e-9: 3.4457, 3200e-9: 3.43438589865999}                                            
    MoveIon=0
    Egap_SI = 3.04*e
    Egap_SI_indirect = 0.98*e
    PulseType  = 3 #STH: 2, StH: 3
    wavelength = 800e-9
    ShowPlot   = 1
    PlotWpiAxis = True
    NumberOfAtoms=2. #Si
    PulseType=3 #StH
    ## Plots effect of the STH Rise Time. We even found an incomplete computations thanks to this. Works super well, then. 
    #RiseCycles = 0 #0: undefined, mixes everything. Basically, we must redefine all runs which had no well defined RiseCycle so far. 
    #Efield     = 1E9
    #plotEffectOfRiseTime(PlottingTable, wavelength, 20e-15, Efield, PulseType, ShowPlot)
    
    ## Effects of CEP on single pulse excitation, current and energy, and harmonic spectrum
    RiseCycles = 2 #NOTE: a number different from 0 must be chosen before publication. 
    #Efield     = 1E9
    #plotEffectOfCEP(PlottingTable, wavelength, 20e-15, Efield, PulseType, RiseCycles, ShowPlot)
    ##plotEffectOfCEP(PlottingTable, wavelength, 10e-15, Efield, PulseType, RiseCycles, ShowPlot)
    #Efield     = 5E9 #NOTE: Saturation may prevent to see any difference. 
    #plotEffectOfCEP(PlottingTable, wavelength, 20e-15, Efield, PulseType, RiseCycles, ShowPlot)
    ##plotEffectOfCEP(PlottingTable, wavelength, 10e-15, Efield, PulseType, RiseCycles, ShowPlot)
    
    #NOTE: === LOCAL FILTERING
    #PlottingTable_filtered = FilterDatabaseRemove(PlottingTable, 30e-15, IndexTable['pulses1'])
    PlottingTable_filtered = PlottingTable
    
    ## Effect of intensity figures
    #plotEffectOfIntensity_CEP(PlottingTable_filtered, wavelength, 20e-15, False, 0.0, PulseType, RiseCycles, ShowPlot)
    #plotEffectOfIntensity_CEP(PlottingTable_filtered, wavelength, 10e-15, False, 0.0, PulseType, RiseCycles, ShowPlot)
    #plotEffectOfIntensity_CEP(PlottingTable_filtered, wavelength, 30e-15, False, 0.0, PulseType, RiseCycles, ShowPlot)
    
    # Compute ionization probability, effect of the pulse duration
    #ShowPlot=1
    

#    for wavelength in [1600e-9]: #, 407.58e-9]: #, 1600e-9, 3200E-9]: #, 1600e-9, 800e-9]:

    # for wavelength in [3200e-9, 1600e-9, 800e-9, 407.58e-9]: #, 1600e-9, 3200E-9]: #, 1600e-9, 800e-9]:
    for wavelength in [1030e-9]:  # , 1600e-9, 3200E-9]: #, 1600e-9, 800e-9]:
      # try:
        plotEffectOfIntensity_PulseDuration(PlottingTable_filtered,  wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI, Egap_SI_indirect, OpticalIndex[wavelength], MoveIon, Functional="TB09")
      # except:
      #  logger.warning(Header+"** Warning: no data could be found, or bug took place for "+str(wavelength*1E9)+" nm.")
    
    #plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI, Egap_SI_indirect, OpticalIndex[wavelength], MoveIon)
    
    #This one construct the ionization probabilities! 
    #wavelength = 483e-9 #Resonant wavelength results for LDA.
    #wavelength = 407.58e-9 #Resonant wavelength results for TB09. 
    #plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, Egap_SI, Egap_SI_indirect) #This one construct the ionization probabilities! 
    #exit()
    #wavelength = 484e-9 #Resonant wavelength results for LDA.
    #plotEffectOfIntensity_PulseDuration(PlottingTable_filtered, wavelength, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot) #This one construct the ionization probabilities! 
    ## READY
    #wavelength = 800e-9
    #plotEffectOfIntensity_wavelengths(PlottingTable_filtered, 20e-15, 0., PlotWpiAxis, wavelength, PulseType, RiseCycles, ShowPlot, Egap_SI) 
    #plotEffectOfIntensity_wavelengths(PlottingTable_filtered, 10e-15, 0., PlotWpiAxis, wavelength, PulseType, RiseCycles, ShowPlot)
    #plotEffectOfIntensity_wavelengths(PlottingTable_filtered, 30e-15, 0., PlotWpiAxis, wavelength, PulseType, RiseCycles, ShowPlot)
    
    ShowPlot=1
    #plotEffectOfWavelength_2dMap(PlottingTable, 20e-15, PulseType, RiseCycles, ShowPlot)
    # plotEffectOfWavelength_Tolerance(PlottingTable, 0.03, 20e-15, 0., PulseType, RiseCycles, ShowPlot, NumberOfAtoms, Egap_SI, Egap_SI_indirect)
    
    #NICE FIGURE
    #plotEffectOfWavelength_Tolerance(PlottingTable_filtered, 1E-4, 20e-15, 0., PulseType, RiseCycles, ShowPlot, NumberOfAtoms, Egap_SI, Egap_SI_indirect, 500e-9, 3300e-9, False, False, True, True, "TB09")
    
    #plotEffectOfWavelength_Tolerance(PlottingTable, 0.03, 30e-15, 0., PulseType, RiseCycles, ShowPlot, Egap_SI)
    ### OLD STUFFS
    #plotEffectOfIntensity_CEP(PlottingTable, wavelength, 10e-15, PulseType, RiseCycles, ShowPlot)
    ## Plots the N_exc(wavelength, field) (2D plot) #NOTE: this can be compared with GrojoPRB2013. 

#}}}

## Prepare all plots for bicolor paper with V. Zhukov [1 Phys. Rev. Lett + 1 Phys. Rev. A]
def plot_Zhukov_BicolorMaps(): #{{{
    Egap_SI = 2.56*e; Egap_SI_indirect=0.51*e
    
    ## Plotting E1, E2 maps 
    for tau in [30e-15]:
        PulseType  = 3; RiseCycles = 2 #StH corrected. 
        #PulseType  = 2;  RiseCycles = 0 #bugged STH. 
        
        ##Bicolor_plotEffectOfCEP(PlottingTable, 1030e-9, tau, 1E9, 515e-9, tau, 1E9, PulseType, RiseCycles)
        ##Bicolor_plotEffectOfCEP(PlottingTable, 1030e-9, tau, 9E9, 515e-9, tau, 9E9, PulseType, RiseCycles)
        ##exit()
        
        ## DONE: 800+1600
        ## TODO: 1030+515, 800+400
        ## TODO: 800+1030
        
        #### These runs are separatly SUB-band gap, but jointly ABOVE band gap. 
        ### Commensurable 800+1600
        #wavelength1 = 800e-9; wavelength2 = 1600e-9
        CEP1=0.; CEP2=0.; Egap=2.56*e; meff=0.2226; 
        ## BicolorE1E2plotDictionnary={'Nexc': 0, 'E2': 4, 'J': 1, 'Energy': 2, 'Wpi': 3}
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['Nexc'])
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['NexcGrad'])
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['J'])
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['Energy'])
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['Wpi'])
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['E2'])
        ### Generate the tables from Zhukov files and plot them
        ##VPZ_Wpi2D(wavelength1, wavelength2, 4E9, 2E9)
        
        ##libKeldyshZhukov.VZ_generateWpiTables(Efield1, Efield2, wavelength1, wavelength2, CEP1, CEP2, Egap, meff)
        ##plot_VZ_Bicolor()
        
        ### Commensurable ## NOTE: these results will never reveal multiphotonic steps, because of free-carrier absorption, but also because of the multiple bands. Only possibility is to make sections for a given field, and to analyze the number of multiphoton absorption. 
        ## BUG: some new runs on Prometheus keep RiseCycle=0 although they should be RiseCycle=2. 
        ## Then, I impose RiseCycles = 0 to avoid selecting them. TODO: fix it when duplicates with various RISE cycles would appear. 
        wavelength1 = 1030e-9; wavelength2 = 515e-9; RiseCycles=0 #NOTE: just to lift selection
        #wavelength1 = 800e-9; wavelength2 = 1030e-9; RiseCycles=0 #NOTE: just to lift selection
        plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['Nexc'])
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['NexcGrad'])
        plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['J'])
        plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['Energy'])
        plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['Wpi'])
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  BicolorE1E2plotDictionnary['E2'])
        ##libKeldyshZhukov.VZ_generateWpiTables(Efield1, Efield2, wavelength1, wavelength2, CEP1, CEP2, Egap, meff)
        ##VPZ_Wpi2D(wavelength1, wavelength2, 0.8E10, 1E10)
        ### Incommensurable
        #wavelength1 = 800e-9 
        #wavelength2 = 1030e-9; 
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  0)
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  1)
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  2)
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  3)
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  4)
        #plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  5)
                
        ### Plotting the Stark effect for 2 colors (1: produces N_exc(omega). 2: Energy(omega).)
        
    # Plotting f(wavelength2) for several fields, and a given wavelength1    
    for tau in [30e-15]: #{{{
        ShowPlot=1
        FieldLimit=1E10 #5E9
        wavelength1 = 1030e-9;
        PulseType  = 3; RiseCycles = 0 #NOTE: RiseCycles=2 would be better #StH corrected. 
        for Field1 in [2.5E8]: #, 1.0E9]: # TODO: how to add the E2=0 & E1!=0 ?
            #try:
            Bicolor_plotEffectOfWavelength_Tolerance(PlottingTable, wavelength1, 1E-5, Field1, FieldLimit, tau, 0., PulseType, RiseCycles, ShowPlot, Egap_SI, Egap_SI_indirect)
            logger.info((Header+"Bicolor_plotEffectOfWavelength_Tolerance() success for (tau, wavelength, Efield1)=("+str(tau*1E15)+" fs, "+str(wavelength1*1E9)+" nm, "+str(Field1*1E-9)+" V/nm)"))
            #except:
                #print(Header+"Bicolor_plotEffectOfWavelength_Tolerance() failed for (tau, wavelength, Efield1)=("+str(tau*1E15)+" fs, "+str(wavelength1*1E9)+" nm, "+str(Field1*1E9)+" V/nm)")
    exit()
    #}}}
    #VPZ_Wpi2D(wavelength1, wavelength2, 0.8E10, 1E10)
    
    ## Stark analysis of the bicolor simulations. 
    # Plot N_exc(E2; E1=fixed)
    
    
    # How to plot absorbed energy as function of the CEP2. 
#}}}

## Prepare all plots for shift current paper with A. Rubio
def plot_ShiftCurrents():
    for tau in [30e-15, 50e-15]:
        PulseType  = 1  #SinSquared
        RiseCycles = 0  #We dont care  
        wavelength1 = 3100e-9
        for wavelength2 in [1.55e-6]:
            ## Commensurable
            plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  0)
            plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  1)
            plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  2)
            plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  3)
            ## Incommensurable                                                                                        
            plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  0)
            plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  1)
            plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  2)
            plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1, wavelength2, tau, PulseType, RiseCycles,  3)
    
def plot_HHG():
        for tau in [20e-15, 30e-15, 50e-15]:
            PulseType = 3 #1: SinSquared. 3: STH
            RiseCycles = 0 #We dont care
            wavelength1 = 3100e-9
            ## Which type of plot for HHG? 
            # HHG(omega; run1, run2, run3, run4, run5) #to check convergency [USEFUL]
            # HHG(omega; LAMBDA1, LAMBDA2, given fields) [RISKY, low visibility]
            # HHG(omega; E1, E2; LAMBDA1, LAMBDA2) [EFFICIENT, VISIBLE, MEANINGFUL] Expected: increase of the piedestal/cutoff
            # HHG(omega; tau1; LAMBDA1) Expected: reduction of the harmonics with time. Gabor transform necessary. 
            # HHG(omega; lambda1) #Figure made by N. Tancogne in his PRL 2016. Great to see the cutoff. 
            # Filter: we don't care about CEP

## Prepares plots for metallic absorption computed via TDDFT            
def plot_Vorobyev_Au(material):
    PlottingTable_filtered = PlottingTable
    #for tau in [30e-15, 60e-15]:
    ShowPlot=1
    RiseCycles=0 #does not matter here
    MoveIon=0
    PlotWpiAxis=False #we are interested in energy
    wavelength1=800e-9 #249e-9 also available
    Egap_SI=0e0; Egap_SI_indirect=0e0 #metal
    OpticalIndex_Mo={1030e-9: 2.4201,  800e-9: 3.6371,  515e-9: 3.7362, 400e-9: 3.030048754875381, 249e-9: 1.49083003952569} #Ordal
    OpticalIndex_Au={1030e-9: 0.24212, 800e-9: 0.15352, 515e-9: 0.71912, 400e-9: 1.657892739273246, 249e-9: 1.48523320158103} #Palik
    OpticalIndex={'Mo': OpticalIndex_Mo, 'Au': OpticalIndex_Au}
    PrimitiveCellArea = {'Au': 6.52E-20, 'Mo': 7.07E-20} #m2
    RatioExternal={'Au': 24.06, 'Mo': 1.91} #origin: In Yabana et al PRB (2012), this E_out / E_in = OpticalIndex**2
    #try:
    if(material=="Au"):
        for PulseType in [libPlotting.PulseTypeID['SinSquared'], libPlotting.PulseTypeID['StH']]:
        #0: Gaussian; 1: SinSquared
            plotEffectOfIntensity_PulseDuration_Metals(PlottingTable_filtered,  wavelength1, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, PrimitiveCellArea[material], OpticalIndex[material][wavelength1], RatioExternal[material], MoveIon)
        #except:
            #print Header+"** Warning: no data could be found, or bug took place for 800 nm."
    elif(material=="Mo"):
        for PulseType in [libPlotting.PulseTypeID['SinSquared']]:
        #0: Gaussian; 1: SinSquared
            plotEffectOfIntensity_PulseDuration_Metals(PlottingTable_filtered,  wavelength1, 0., PulseType, RiseCycles, PlotWpiAxis, ShowPlot, PrimitiveCellArea[material], OpticalIndex[material][wavelength1], RatioExternal[material], MoveIon)
    else:
        logger.error((Header+"** Error: material is not well specified."))
        exit()

