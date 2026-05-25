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

## @package plotFinalQuantities
# Backward-compatible entry point.  All functions have been moved to focused
# sub-modules; this file re-exports them so that existing import statements
# continue to work unchanged.
#
# Sub-modules:
#   _plotCommon          -- shared state & database loading
#   plotSinglePulse      -- single-colour / single-pulse plots
#   plotWavelength       -- wavelength-dependent plots
#   plotDoublePulse      -- two-colour / double-pulse plots
#   plotPaperCollections -- paper-specific figure assemblers

from octopus_slabs.Libs._plotCommon import *          # noqa: F401,F403
from octopus_slabs.Libs.plotSinglePulse import *      # noqa: F401,F403
from octopus_slabs.Libs.plotWavelength import *       # noqa: F401,F403
from octopus_slabs.Libs.plotDoublePulse import *      # noqa: F401,F403
from octopus_slabs.Libs.plotPaperCollections import * # noqa: F401,F403

# VALID but some degeneracies !
#plotEffectOfIntensity(PlottingTable, 10e-15, 1, 1, 1) #validated
#plotEffectOfIntensity_wavelengths(PlottingTable, 20e-15, 1, 0) #validated, degeneracies are due to undefined RiseTime
#plotEffectOfIntensity_wavelengths(PlottingTable, 30e-15, 0) #validated
#plotEffectOfIntensity_wavelengths(PlottingTable, 40e-15, 0) #validated
#plotEffectOfIntensity_wavelengths(PlottingTable, 50e-15, 0) #validated
#plotEffectOfIntensity_PulseDuration(PlottingTable, 800e-9, 1) #TODO: remove the different CEPs and keep only tau dependence. 

# TYPICAL RUNS
#plotEffectOfWavelength_2dMap(PlottingTable) #VALID works, but not enough data for now. 
#plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple
#plotEffectOfDoublePulse_Color1_Color2_SingleFieldStrength(PlottingTable, 5E9, 0) #VALID, but don't apply to TOO LARGE wavelength sets
#plotEffectOfDoublePulse_Color1_Color2_SingleFieldStrength(PlottingTable, 5E9, 1) #VALID, but don't apply to TOO LARGE wavelength sets
#plotEffectOfDoublePulse_Color1_Color2_SingleFieldStrength(PlottingTable, 5E9, 2) #VALID, but don't apply to TOO LARGE wavelength sets
#plotEffectOfDoublePulse_Color1_Color2_SingleFieldStrength(PlottingTable, 5E9, 3) #VALID, but don't apply to TOO LARGE wavelength sets

# DEVELOPMENT RUNS
#plotEffectOfWavelength_Tolerance(PlottingTable) #considers a tolerance to plot various fields
#plotEffectOfCEP(PlottingTable)
# plotEffectOfFrozenXC(PlottingTable)

#plotEffectOfDoublePulse_Efield1_Efield2_batchAllColors(PlottingTable) #BUGGY
#plotEffectOfDoublePulse_sameAmplitude_samePulseDurations(PlottingTable)
#plotEffectOfDoublePulse_diffWavelengths(PlottingTable)

# Execute only if called directly from console
if(__name__=="__main__"):
    #====================================== EXECUTION ZONE ============================

    print("Welcome to Octopus-slab plotter.")
    print("================================\n")
    print("** Info: Automatic detection of the routine to run...")

    dirname=os.path.realpath(".") #./02-Slab2D")  ##TODO: this is no longer valid. 

    logger.info("** Info: directory name: "+str(dirname))

    ## Zhukov single color plotting [Phys. Rev. B]

    if "BicolorSummaries" in dirname:
        ## Zhukov bicolor plotting
        logger.info(Header+"** Info: running bicolor analyses of Si LDA...")
        plot_Zhukov_BicolorMaps() #LDA gap
        
    elif "03-Bulk3D" in dirname:
        ## NOTE: this is rather for forcing usage
        #plot_Zhukov_SingleColor_LDA()
        # plot_Zhukov_SingleColor_TB09()
        ## NOTE: this is fine for local usage
    #    exit()
        # if "TB09" in dirname: #WARNING: TODO: This criterion is no longer applicable,
        # since database construction now contains functionals. Better solution would be to detect the
        # available functional, to map the gap and effective masses to a database.
        # Is NOMAD capable to do that? And ASE?
        logger.info("Running LDA single-color Keldysh analyses...")
        plot_Zhukov_SingleColor_LDA()
        # logger.info("Detected TB09 computations for Si...")
        # plot_Zhukov_SingleColor_TB09()

    elif "09-Au" in dirname: 
        logger.info((Header+"** Info: detected computations for Au..."))
        plot_Vorobyev_Au('Au')
    elif "06-Mo" in dirname: 
        logger.info((Header+"** Info: detected computations for Mo..."))
        plot_Vorobyev_Au('Mo')
    elif "14-SiO2" in dirname: 
        logger.info((Header+"** Info: detected simulation for SiO2..."))
        plot_Zhukov_SingleColor_SiO2()
    else: 
        logger.warning("Other routines are available, but may not be up to date.")
        
    ## Rubio shift current plotting
    #plot_ShiftCurrents()

    ## Preparing HHG spectra dependencies for Hamed Merdji
    #plot_HHG()

    ## Preparing plots for Vorobyev

    logger.info("End of program.")
    #exit()

