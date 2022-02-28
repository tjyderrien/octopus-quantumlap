# if executed from terminal, the code checks which directory is being used and uses parameters for particular material accordingly

from pathlib import Path
import sys
sys.path.append('/home/kristyna.gazdova/octopus-slabs/octopus_slabs/Libs')

# Other options:
# sys.path.insert(1, '/home/kristyna.gazdova/octopus-slabs/octopus_slabs/Libs')
# add the path to PYTHONPATH environment variable

import plotFloquetBandstructure.py as floquet

## PROGRAM BEHAVIOR
SingleKpoint = 0  # Selector for direct execution. 
# 0: plotting the whole band structure for a certain field. 
# 1: plot one single kpoint for many fields and wavelengths.

field_polarization_dir = 0  # 0: x, 1: y, 2: z

## CHOICE OF OUTPUTS
PlotMatrixElements = True
PlotFloquetBandStructure = True  # Provide field dependent band structure at Gamma Point
PlotRabiFrequencyMatrices = False  # for studying Wannier-Stark localization, i.e. suppression of tunneling with zeros of the Bessel function
ShowPlot = False
StorePlot = True #we generate so many figures with this... Let's avoid this at the moment. #starts to be too consuming (800 MB per EPS file...)

Enable_A2 = False #take A2 term into Floquet hamiltonian. NOTE: in ipython notebook, this is not relevant. From console, it is mandatory.

if(__name__=="__main__"):  #runs only if file is called directly
    directory = Path.cwd()
    logger.info("** Info: the current working directory is: ", directory.name) 
    
    if "Au" in directory.name:
        # specific parameters for Au
        num_wavelengths = 500
        wavelength_SI_min = 400e-9
        wavelength_SI_max = 3200e-9
        #SelectedWavelengths = [2000e-9]

        num_fields = 300  # number of electric field values for the plotting
        Efield_SI_min = 1E8 # IntensityToField(1E10 * 1E4)
        Efield_SI_max = 1E12 # IntensityToField(1E19 * 1E4)
        #SelectedFieldsForPlotting= [0.001E9, 0.01E9, 0.05E9]  #selectors dont overlap with num_fields (both are used, though! You'll see :) ).

        num_time_steps = 300 #should be small enough to ensure convergence

        # Warning: See Yabana PRB 2012 + Abraham-Minkowski momentum conservation
        UseRefractiveIndex = True  #True: limits the plots maximum as function of wavelength (makes field wavelength-dependent)
        PlotIntensity = True  # True: plot with I. False: plot with E.

        SizeOfDots = 0.5  # purely cosmetic!
        SizeX = 8.2; SizeY=6.

        LogPlotting = 1 #0: regular linear scale, 1: regular logscale
        verbose = 0  # Show info and debug messages.

        BesselRootsID_max = 10
        
        logger.info("** Plotting Floquet Band structure for Au")
        
        if (SingleKpoint == 1):  # for now, SingleKpoint=0 in floquet -> if you wish to plot SingleKpointManyFields, set SingleKpoint to 1 in floquet
            floquet.SingleKpointManyFields(Enable_A2)
        else:
            for wavelength in [2000e-9]: #1030e9 #SelectedWavelengths:
                # for field in [1E6, 1E7, 5E7]:
                for field in [0.05E9]: #[1E9, 2E9, 3E9, 4E9, 5E9, 6E9, 7E9, 8E9, 9E9]: #[1E9, 2E9, 3E9, 4E9, 5E9]
                    floquet.MultipleKpointsOneField(field, wavelength, Enable_A2) 
    elif "Si" in directory.name:
        # specific parameters for Si
        num_wavelengths = 500
        wavelength_SI_min = 400e-9
        wavelength_SI_max = 3200e-9
        #SelectedWavelengths = [2000e-9]

        num_fields = 300  # number of electric field values for the plotting
        Efield_SI_min = 1E8 # IntensityToField(1E10 * 1E4)
        Efield_SI_max = 1E12 # IntensityToField(1E19 * 1E4)
        #SelectedFieldsForPlotting= [0.001E9, 0.01E9, 0.05E9]  #selectors dont overlap with num_fields (both are used, though! You'll see :) ).

        num_time_steps = 300 #should be small enough to ensure convergence

        # Warning: See Yabana PRB 2012 + Abraham-Minkowski momentum conservation
        UseRefractiveIndex = True  #True: limits the plots maximum as function of wavelength (makes field wavelength-dependent)
        PlotIntensity = True  # True: plot with I. False: plot with E.

        SizeOfDots = 0.5  # purely cosmetic!
        SizeX = 8.2; SizeY=6.

        LogPlotting = 1 #0: regular linear scale, 1: regular logscale
        verbose = 0  # Show info and debug messages.

        BesselRootsID_max = 10
        
        logger.info("** Plotting Floquet Band structure for Si")
        
        if (SingleKpoint == 1):  # for now, SingleKpoint=0 in floquet -> if you wish to plot SingleKpointManyFields, set SingleKpoint to 1
            floquet.SingleKpointManyFields(Enable_A2)
        else:
            for wavelength in [2000e-9]: #1030e9 #SelectedWavelengths:
                # for field in [1E6, 1E7, 5E7]:
                for field in [0.05E9]: #[1E9, 2E9, 3E9, 4E9, 5E9, 6E9, 7E9, 8E9, 9E9]: #[1E9, 2E9, 3E9, 4E9, 5E9]
                    floquet.MultipleKpointsOneField(field, wavelength, Enable_A2)
    elif "SiO2" in directory.name:
        # specific parameters for SiO2
        num_wavelengths = 500
        wavelength_SI_min = 400e-9
        wavelength_SI_max = 3200e-9
        #SelectedWavelengths = [2000e-9]

        num_fields = 300  # number of electric field values for the plotting
        Efield_SI_min = 1E8 # IntensityToField(1E10 * 1E4)
        Efield_SI_max = 1E12 # IntensityToField(1E19 * 1E4)
        #SelectedFieldsForPlotting= [0.001E9, 0.01E9, 0.05E9]  #selectors dont overlap with num_fields (both are used, though! You'll see :) ).

        num_time_steps = 300 #should be small enough to ensure convergence

        # Warning: See Yabana PRB 2012 + Abraham-Minkowski momentum conservation
        UseRefractiveIndex = True  #True: limits the plots maximum as function of wavelength (makes field wavelength-dependent)
        PlotIntensity = True  # True: plot with I. False: plot with E.

        SizeOfDots = 0.5  # purely cosmetic!
        SizeX = 8.2; SizeY=6.

        LogPlotting = 1 #0: regular linear scale, 1: regular logscale
        verbose = 0  # Show info and debug messages.

        BesselRootsID_max = 10
        
        logger.info("** Plotting Floquet Band structure for SiO2")
        
        if (SingleKpoint == 1):  # for now, SingleKpoint=0 in floquet -> if you wish to plot SingleKpointManyFields, set SingleKpoint to 1 in floquet
            floquet.SingleKpointManyFields(Enable_A2)
        else:
            for wavelength in [2000e-9]: #1030e9 #SelectedWavelengths:
                # for field in [1E6, 1E7, 5E7]:
                for field in [0.05E9]: #[1E9, 2E9, 3E9, 4E9, 5E9, 6E9, 7E9, 8E9, 9E9]: #[1E9, 2E9, 3E9, 4E9, 5E9]
                    floquet.MultipleKpointsOneField(field, wavelength, Enable_A2)