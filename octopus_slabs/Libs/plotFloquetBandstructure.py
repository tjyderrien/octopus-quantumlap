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
# along with this program. If not, see <http://www.gnu.org/licenses/>

## @package plotFloquetBandstructure
# Provides routines and prepare dressed band structure from ground state matrix elements of a solid crystal.
# Warning: symmetry breaking in taken in effective manner, and not explictly in dipolar matrix elements.
"""
Explanation of the number of method employed by Octopus for indexing the k-grid files:
- Numbering of the files is N_files = (kgrid) + (kpoint path)
- 1,2,3 directions for dipolar matrix elements (x,y,z), added a suffix of the dipolar files ".x###".
- Example: k-grid size=(16 x 16 x 16) = 4096 files for 3D k-band structure.
- One k-point path is given at the end of the numbering. e.g. 1001 files remain for k-point path.
"""

#import pdb
from matplotlib.cm   import get_cmap
from cycler import cycler
from scipy.constants import c, epsilon_0, h, e
from scipy.special   import jn_zeros as BesselJZeros
import os
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

import octopus_slabs.Libs.libPlotting as libPlotting
import octopus_slabs.Libs.libString   as libString
import octopus_slabs.Libs.libAtomicUnits     as au
from   octopus_slabs.Libs.libLogging import init_logger

from spp_extended_theory.Keldysh.libKeldyshZhukov import IntensityToField
from spp_extended_theory.Keldysh.libStark import ComputeFloquetBandStructure
import spp_extended_theory.Libs.libImportOpticalData as libOpticalData

## PROGRAM BEHAVIOR
SingleKpoint = 1  # Selector for direct execution. 
# 0: plotting the whole band structure for a certain field. 
# 1: plot one single kpoint for many fields and wavelengths.

field_polarization_dir = 0  # 0: x, 1: y, 2: z

## CHOICE OF OUTPUTS
PlotMatrixElements = True
PlotFloquetBandStructure = True  # Provide field dependent band structure at Gamma Point
PlotRabiFrequencyMatrices = False  # for studying Wannier-Stark localization, i.e. suppression of tunneling with zeros of the Bessel function
ShowPlot = False
StorePlot = True #we generate so many figures with this... Let's avoid this at the moment. #starts to be too consuming (800 MB per EPS file...)
Enable_A2 = False #take A2 term into Floquet hamiltonian

## CHOOSING INPUTS

num_wavelengths = 500
wavelength_SI_min = 400e-9
wavelength_SI_max = 3200e-9
SelectedWavelengths = [400e-9, 515e-9, 800e-9, 1030e-9, 3200e-9]

num_fields = 200  # number of electric field values for the plotting
Efield_SI_min = 1E9 # IntensityToField(1E10 * 1E4)
Efield_SI_max = 5E11 # IntensityToField(1E19 * 1E4)
SelectedFieldsForPlotting=[0.5E9, 1E9, 5E9, 1E10] #selectors dont overlap with num_fields (both are used, though! You'll see :) ).

num_time_steps = 300 #should be small enough to ensure convergence

UseRefractiveIndex = False  # True: limits the plots maximum as function of wavelength (makes field wavelength-dependent)
PlotIntensity = True  # True: plot with I. False: plot with E.

SizeOfDots = 0.5  # purely cosmetic!
SizeX = 8.2; SizeY=6.

LogPlotting = 1 #0: regular linear scale, 1: regular logscale
verbose = 0  # Show info and debug messages.

BesselRootsID_max = 10

########### DONT TOUCH BELOW, thanks! #############
logger = init_logger(__name__, verbose=False) #"plotFinalQuantities")

## Building list of files

# Filenames = "ks_me_multipoles.*"
Filenames = "static/ks_me_dipole.*" #name of dipolar matrix elements in octopus
bandstructure = "bandstructure"  #name of bandstructure file in octopus NOTE: file may be absent sometimes
LocalFolder = "./"
filename = Filenames.split('.')[0]
# print Filenames
# 1. Capture the static/dipole files.
listoffiles = libString.GetOrderedListOfFiles(Filenames)
del Filenames
# print listoffiles


## Building materials dictionnary
dirname=os.path.realpath(".") #catch the present directory absolute path

logger.info(("** Info: directory name: "+str(dirname)))

logger.info("** Info: capturing the present material...")
try:
    if "03-Bulk3D" in dirname:
        material = "Si"
    elif "09-Au" in dirname:
        material = "Au"
    elif "06-Mo" in dirname:
        material = "Mo"
    elif "14-SiO2" in dirname:
        material = "SiO2"
    logger.info("** Info: detected the material "+material)
except:
    logger.error(Header+"Which material are you trying to analyze?")
    exit()


spin_occupation = 2  # 1: levels with same spin should split, 2: levels of different spins don't split.
Complex_valued_ME = False #True for non-Hermitian eigenvalues?

Efield_AU_min = au.Field_SI_to_AU(Efield_SI_min)
Efield_AU_max = au.Field_SI_to_AU(Efield_SI_max)

OpticalIndex_Si = {
    # NOTE: these are experimental data for Si (from Palik). They should be from linear TDDFT response !
    343e-9: 5.2749, 400e-9: 5.5674, 483e-9: 4.4028, 484e-9: 4.3964,
    800e-9: 3.6924167231905, 1600e-9: 3.4826,
    2200e-9: 3.4548, 2600e-9: 3.4457, 3200e-9: 3.43438589865999, 10E-6: 3.4181}

OpticalIndex_SiO2 = {
    # NOTE: these are experimental data for SiO2 (from refractiveindex.info). They should be from linear TDDFT response !
    343e-9: 1.5669, 400e-9: 1.5577, 483e-9: 1.5499, 484e-9: 1.5498,
    800e-9: 1.5383, 1600e-9: 1.5270,
    2200e-9: 1.5175, 2600e-9: 1.5094, 3200e-9: 1.4996}

OpticalIndex_Mo = {3200e-9: 1.6626, 1600e-9: 1.5686, 1030e-9: 2.4201, 800e-9: 3.6371, 400e-9: 0.42388}

OpticalIndex_Au = {1600e-9: 0.55401, 1030e-9: 0.24212, 800e-9: 0.15352, 400e-9: 1.4684, 343e-9: 1.4809}

# In band gap materials, we are interested in plotting evolution of energy levels at Gamma point as function of laser field.
Si = {'MPI_numbers': [1], 'kpoints_manual': 16,
      'nb_atoms': 1, 'unocc_states': 4,
      'Z_electrons': 4, 'OpticalIndex': OpticalIndex_Si, 'OpticalIndexFile': "Si-Palik",
      'xmin_xmax': (2E11, 7E16),
      'ymin_ymax': (-5, 15), 'CaptureDirectGap': True, 'KpointPath_len': 301,
      'K_labels': [r'$\Gamma$', r'$X$', r'$U|K$', r'$\Gamma$', r'$L$', r'$W$', r'$X$']}

SiO2 = {'MPI_numbers': [1], 'kpoints_manual': 8,
        'nb_atoms': 1, 'unocc_states': 6,
        'Z_electrons': 24, 'OpticalIndex': OpticalIndex_SiO2, 'OpticalIndexFile': "SiO2-Palik",
        'xmin_xmax': (2E10, 7E16),
        'ymin_ymax': (-8, 16), 'CaptureDirectGap': True, 'KpointPath_len': 301,
        'K_labels': [r'$\Gamma$', r'$M$', r'$K$', r'$\Gamma$', r'$A$', r'$L$', r'$H$', r'$A|L$', r'$M|H$', r'$K$']}

# In metals, we are interested in plotting band structure on the full path of high symmetry k-points
Mo = {'MPI_numbers': [1], 'kpoints_manual': 16,
      'nb_atoms': 1, 'unocc_states': 3,
      'Z_electrons': 8, 'OpticalIndex': OpticalIndex_Mo, 'OpticalIndexFile': "Mo-Palik",
      'xmin_xmax': (2E10, 7E16),
      'ymin_ymax': (-30, 30), 'CaptureDirectGap': False, 'KpointPath_len': 601,
      'K_labels': [r'$\Gamma$', r'$H$', r'$N$', r'$\Gamma$', r'$P$', r'$H|P$', r'$N$'], 
      'Fermi_Energy': 12.28
}

Au = {'MPI_numbers': [1], 'kpoints_manual': 16,
      'nb_atoms': 1, 'unocc_states': 3,
      'Z_electrons': 14 - 3, 'OpticalIndex': OpticalIndex_Au, 'OpticalIndexFile': "Au-Palik",
      'ymin_ymax': (-2, 35), 'CaptureDirectGap': False, 'KpointPath_len': 1001,
      'xmin_xmax': (2E10, 7E16),
      'K_labels': [r'$\Gamma$', r'$X$', r'$U|K$', r'$\Gamma$', r'$L$', r'$W$', r'$X$'], 
      'Fermi_Energy': 8.58 #eV at 300K
}
      #'K_labels': [r'$\Gamma$', r'$X$', r'$W$', r'$U|K$', r'$\Gamma$', r'$L$', r'$U$', r'$W$', r'$L$', r'$K|U$', r'$X$']

Materials = {'Si': Si, 'SiO2': SiO2, 'Mo': Mo, 'Au': Au}  # this contains all materials properties we need to automatize

# 2. Just shortening the filenames here
nb_atoms = Materials[material]['nb_atoms']  # 1 #number of atoms #we actually treat a bulk, so it's one.
Z_electrons = Materials[material][
    'Z_electrons']  # SiO2: 24 electrons. #Si: 4 #number of modeled valence electrons in the simulation box
# (3p)² x 2 for Si (biatomic) lattice.
# SiO2 ? matrix is 30 e- states, 6 unocc.
unocc_states = Materials[material][
    'unocc_states']  # 4 #SiO2: 6. #Si: 4. This is defined in the inp. If not, then it's zero.
OpticalIndexFile = Materials[material]['OpticalIndexFile'] #routes to continous optical data from Palik

CaptureDirectGap = Materials[material]['CaptureDirectGap']

MPI_numbers = Materials[material]['MPI_numbers']  # [1, 2, 3]
kpoints_manual = Materials[material]['kpoints_manual']  # 16 #number of kpoints in the grid
# OpticalIndex = Materials[material]['OpticalIndex']
KpointPath_len = Materials[material]['KpointPath_len']
# print
K_labels = Materials[material]['K_labels']

# photon_energy_SI_t=np.linspace(h*c/wavelength_SI_max/e,h*c/wavelength_SI_min/e, num_wavelengths) #eV
wavelength_SI_t = np.linspace(wavelength_SI_min, wavelength_SI_max, num_wavelengths)
# wavelength_SI_t = h*c / e / photon_energy_SI_t #meters
# wavelength_SI_t = wavelength_SI_t[::-1] #to keep it in increasing order

logger.info(wavelength_SI_t)

## Prepare energy diagram as function of laser field strength at a single K point #typical for M. Stockman papers
# @param k_point_index_of_interest: 1 is $\Gamma$. For other k-points, think !
def SingleKpointManyFields(k_point_index_of_interest=1, wavelength_t=wavelength_SI_t):  # {{{
    Header="plotFloquetBandStructure.py: SingleKpointManyFields(): "
    global Efield_SI_min, Efield_SI_max, kpoints_manual, Z_electrons, unocc_states, OpticalIndexFile, Enable_A2, \
        SelectedWavelengths, SelectedFieldsForPlotting

    plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y'])))

    ## Capturing the Gamma k-point file in results of Octopus
    # TODO: generalize to a desired k-point. For now it is using GAMMA point. Maybe link to ASE can do that? Or a lib of symmetries.
    try:
        kpoints_number_inp = libString.LookUpStringInFile("info",
                                                          "Dimensions of the k-point grid")  # WARNING: it seems this has changed.
        for i in np.arange(0, 10):  # TODO: replace with .strip() ?
            kpoints_number_inp = kpoints_number_inp.replace("  ", " ")
            if (verbose == 1):
                logger.info("** Info: decreasing number of spaces")
                logger.info(kpoints_number_inp)
        kpoints_number_array = kpoints_number_inp.split("=")
        kpoints_number_array = kpoints_number_array[1].split(" ")
        kpoints = int(kpoints_number_array[1])
        logger.info(kpoints)
    except IndexError: #recent octopus versions removed some info from info file
        kpoints = kpoints_manual  # manual input of number of k-grid points in 1 direction (assumes same resolution in 3 directions)
        logger.error("** Error: failed to capture the number of k points automatically. Manual input is employed: kpoints=" + str(
            kpoints))
        # exit()

    logger.info("Info: Capturing the file for the transition of interest.")
    #k_point_index_of_interest=1 #1 is GAMMA point. Points are sorted by distance to Gamma point. 
    # WARNING: this captures the GAMMA POINT. 
    gamma_point_files, gamma_point_files_index = libString.ExtractFileForGivenKpoints(filename, kpoints, listoffiles, [k_point_index_of_interest],
                                                                                      CaptureDirectGap, verbose,
                                                                                      field_polarization_dir)
    #breakpoint()

    matrixelements, matrix_side = libString.ReadDipolarMatrixElements(Z_electrons, unocc_states, gamma_point_files,
                                                                      Complex_valued_ME)

    ## 5. Extract eigen values at k-point <kpoints_index> from "info" file #WARNING: selecting first k-point (Gamma) here!
    eigenvalues = libString.GetEigenValue(gamma_point_files_index[0], Z_electrons + unocc_states, "info", verbose)
    if (verbose == 1):
        logger.debug(eigenvalues)

    # tunneling_stops=[]

    epsilon_Si0       = libOpticalData.ImportPalikDatabase_epsilon_fromNK(SelectedWavelengths, OpticalIndexFile)  ## WARNING: replacing by LDA is DIFFICULT, since PRL Rubio is defined down to 2 eV. Nothing below!! #TODO: this should be material-dependent. Indicate path to fetch the right file in library at the top.
    RefractiveIndex_t = np.real(np.sqrt(epsilon_Si0)) #This is a continuum formulation from Palik data
    if (LogPlotting==1):
        Efield_SI_log = np.linspace(np.log10(Efield_SI_min), np.log10(Efield_SI_max), num_fields)
        Efield_SI = np.power(10., Efield_SI_log)
    elif (LogPlotting==0):
        Efield_SI = np.linspace(Efield_SI_min, Efield_SI_max, num_fields)
    else: 
        Efield_SI = SelectedFieldsForPlotting

    SizeOfReplica = (Z_electrons+unocc_states)*( 2 * MPI_numbers[0] + 1 )
    eigen_idx = np.arange(0, SizeOfReplica ) #we only care about simple replicas at the moment
    # Creating data array for easier visualization
    replicas_data = xr.DataArray(np.zeros((len(Efield_SI), len(SelectedWavelengths), len(eigen_idx))),
                                 dims=("E", "wavelength", "eigen_idx"),
                                 coords={"E": Efield_SI, "wavelength": SelectedWavelengths, "eigen_idx": eigen_idx})
    # TODO: replicas_data should be revised completely. Fig1: contains few wavelengths, many fields. Fig2: many
    #  wavelengths, few fields. replicas_data should have different dimensions for each.
    logger.info("Prepared tensor for storing the results. ")
    if(not PlotFloquetBandStructure): 
        logger.warning("If you want to visualize Floquet band structures, check the value of PlotFloquetBandStructure!")
    if (PlotFloquetBandStructure):
        # for wavelength_SI in wavelength_SI_t:
        
        # replica_t = np.zeros(( num_wavelengths, num_fields )) #tenseur second order
        for wavelength_idx in np.arange(0,len(SelectedWavelengths)): #altered from original to get figures for desired wavelengths
            #wavelength_SI = wavelength_t[wavelength_idx] #original
            #wavelength_t[wavelength_idx] = SelectedWavelengths #one of many unsuccessful attempts to make it work
            wavelength_SI = SelectedWavelengths [wavelength_idx]
            OpticalIndex = RefractiveIndex_t [wavelength_idx]    #This is a continuum formulation from Palik data


            ## Renormalize the intensity range using refractive index - TODO: re-enable using xarray
#            if (UseRefractiveIndex):
#                Efield_SI_max = Efield_SI_max / np.sqrt(OpticalIndex) #field becomes wavelength dependent?

            ## Plot the splitted energy levels at Gamma.
            for MPI_number in MPI_numbers:
                fig1, ax1 = plt.subplots()
                # for field in Efield_SI:
                for field_idx in np.arange(0,len(Efield_SI)):
                    field = Efield_SI[field_idx]
                    eigenvalues, replicas, dipoles = ComputeFloquetBandStructure(filename, nb_atoms, Z_electrons,
                                                                                 kpoints, unocc_states, MPI_number,
                                                                                 field, wavelength_SI, num_time_steps,
                                                                                 eigenvalues, matrixelements,
                                                                                 Complex_valued_ME,
                                                                                 field_polarization_dir,
                                                                                 spin_occupation, verbose, Enable_A2)
                    # Plotting using external intensity, whence internal intensity acts on the sample
                    if (PlotIntensity):
                        unitx = 1E-4
                        field = 0.5 * c * epsilon_0 * field ** 2 * OpticalIndex
                    else:
                        unitx = 1E0
                        field = field * np.sqrt(OpticalIndex)
                    ax1.scatter(np.ones(np.size(replicas)) * field * unitx, au.Energy_Hartree_to_eV(np.real(replicas)),
                                c="red", s=SizeOfDots,
                                label="Floquet replicates")  # these are the interesting lines for the final plot
                    # ax2.scatter(np.ones(np.size(replicas)) * wavelength_SI * 1E9, au.Energy_Hartree_to_eV(np.real(replicas)), c="red", s=SizeOfDots, label="Floquet replicates")

                    # if(MPI_number==1):
                    replicas_data[field_idx, wavelength_idx, :] = np.real(replicas) #array of lists, we might need XARRAY
                        # breakpoint()
                    # TODO: Prepare the table [wavelength, n, E, Rabi Freq such as BesselJ(n,RabiFreqs_SI/omegaLaser)=0. #BruteForce method. 
                    # TODO: Do we need the function BesselJZeros(n,x)? To search omegaLaser that induce a stop in tunneling! 
                    # tunneling_stops.append(np.sqrt(dipoles * field_AU / BesselJZeros(n,k)))

                ## Plotting eigenvalues ( intensity )
                if (PlotIntensity):
                    Intensity_min = 0.5 * c * epsilon_0 * Efield_SI[0] ** 2 * OpticalIndex
                else:
                    Intensity_min = Efield_SI[0] * np.sqrt(OpticalIndex)

                # Plotting eigenvalues at rest
                ax1.scatter(np.ones(np.size(eigenvalues)) * Intensity_min * unitx,
                            au.Energy_Hartree_to_eV(np.real(eigenvalues)), c="black", s=10, label="GS eigenvalues")

                if (LogPlotting):
                    ax1.set_xscale('log')

                ax1.set_title(r"$n=$" + str(round(MPI_number, 0)) + r", $\lambda=$" + str(
                    int(wavelength_SI * 1E9)) + " nm, k=" + str(kpoints) + r"$^3$")
                figurename = "FloquetStarkShift-" + str(round(wavelength_SI * 1E9, 0)) + "nm-k" + str(
                    kpoints) + "MPI-" + str(round(MPI_number, 0))
                if (not LogPlotting):
                    figurename += "NoLog"
                if (not PlotIntensity):
                    figurename += "-Field"
                    ax1.set_xlabel("Field (V/m)")
                else:
                    ax1.set_xlabel(r"Laser intensity (W/cm$^{2}$)")
                
                if Enable_A2:
                    logger.warning("A2 enabled: ")
                    figurename+= "-EnabledA2"
                ax1.set_ylabel(r"Energy levels (eV)")
                ax1.set_ylim(Materials[material]['ymin_ymax'])
                ax1.set_xlim(Materials[material]['xmin_xmax'])
                try:
                    ax1.axhline(y=Materials[material]['Fermi_Energy'], color='#b3b3b3', linestyle=':') #adds Fermi energy (should be limited to metals: we plot only if the entry exists in dictionnary)
                except KeyError:
                    logger.error("Fermi level will not be plotted. If you are working with a dielectric, there is no drama.")
                # plt.ylim((-8,8)) #SiO2
                fig1.tight_layout()
                if(StorePlot):
                    fig1.savefig(figurename + ".eps")
                    fig1.savefig(figurename + ".png")
                # Plotting dressed levels directly on the figure, from n=1 to MPI_number
                for MPI in np.arange(1, MPI_number + 1):
                    ax1.scatter(np.ones(np.size(eigenvalues)) * Intensity_min * unitx,
                                au.Energy_Hartree_to_eV(np.real(eigenvalues)) + MPI * h * c / wavelength_SI / e, c="red",
                                s=10, label="non-interacting replicates")
                    ax1.scatter(np.ones(np.size(eigenvalues)) * Intensity_min * unitx,
                                au.Energy_Hartree_to_eV(np.real(eigenvalues)) - MPI * h * c / wavelength_SI / e, c="blue",
                                s=10, label="non-interacting replicates")
                if(StorePlot):
                    fig1.savefig(figurename + "-Suppl" + ".eps")
                    fig1.savefig(figurename + "-Suppl" + ".png")
                    logger.info("Output was stored in EPS and PNG files.")
                plt.close(fig1) #free the memory?
                del fig1, ax1
                # plt.legend()
            # plt.show()
            ## Plotting eigenvalues of wavelength
        # End of wavelength loop
        del replicas_data
        # Prepare the plot as function of wavelength for few intensities
        # ax2.set_xlabel("Pump wavelength (nm)")
        

        # ax2.set_xscale("log")
        # We need to reinitialize the table for fig2.
        replicas_data = xr.DataArray(np.zeros((len(SelectedFieldsForPlotting), len(wavelength_t), len(eigen_idx))),
                                     dims=("E", "wavelength", "eigen_idx"),
                                     coords={"E": SelectedFieldsForPlotting, "wavelength": wavelength_t,
                                             "eigen_idx": eigen_idx})
        # replicas_data(Efield_SI, wavelength_t)
        # Select eigenvalues where field is above a certain value
        # breakpoint()

        ## Preparing interpolation for optical data as function of wavelength
        ## WARNING: replacing by LDA is DIFFICULT, since PRL Rubio is defined down to 2 eV. Nothing below!! #TODO: this should be material-dependent. Indicate path to fetch the right file in library at the top.
        epsilon_Si0 = libOpticalData.ImportPalikDatabase_epsilon_fromNK(wavelength_t, OpticalIndexFile)
        RefractiveIndex_t = np.real(np.sqrt(epsilon_Si0))  # This is a continuum formulation from Palik data

        """# Rotation on colors
        cm = get_cmap('hsv')
        NumberOfElements = len(SelectedFieldsForPlotting)
        colors = [cm(1. * i / NumberOfElements) for i in range(NumberOfElements)]
        xy = list(range(NumberOfElements))
        colorlist = [colors[x] for x in xy]
        """
        
        color_idx=0
        for field_idx in np.arange(0,len(SelectedFieldsForPlotting)):
            field = SelectedFieldsForPlotting[field_idx]
            for wavelength_idx in np.arange(0,len(wavelength_t)):
                wavelength_SI = wavelength_t [wavelength_idx]
                OpticalIndex = RefractiveIndex_t [wavelength_idx]

                eigenvalues, replicas, dipoles = ComputeFloquetBandStructure(filename, nb_atoms, Z_electrons,
                                                                                     kpoints, unocc_states, MPI_number,
                                                                                     field, wavelength_SI, num_time_steps,
                                                                                     eigenvalues, matrixelements,
                                                                                     Complex_valued_ME,
                                                                                     field_polarization_dir,
                                                                                     spin_occupation, verbose, Enable_A2)
        
                replicas_data[field_idx, wavelength_idx, :] = np.real(replicas)

            ## One Fig2 per field in SelectedFields
            fig2, ax2 = plt.subplots() #this will contain the wavelength dependence
            ax2.set_prop_cycle(cycler('color', ['c', 'm', 'y', 'k']))
            FloquetKeldysh = replicas_data.sel(E=field, method="nearest") #not essential but good practice for future
            FloquetKeldysh_X = np.repeat(FloquetKeldysh.wavelength, SizeOfReplica)
            FloquetKeldysh_X = h*c/FloquetKeldysh_X / e #eV

            FloquetKeldysh_Y = np.reshape(au.Energy_Hartree_to_eV(FloquetKeldysh.values), (1,np.size(FloquetKeldysh_X)))
            logger.debug(np.shape(FloquetKeldysh_X), np.shape(FloquetKeldysh_Y))
            ax2.scatter(FloquetKeldysh_X, FloquetKeldysh_Y[0], s=SizeOfDots,
                                    label="Floquet replicates")
            ax2.set_xlabel("Pump photon energy (eV)")
            ax2.set_ylabel("Energy levels (eV)")
            ax2.set_title("Field "+str(round(field,1)*1E-9)+" V/nm")
            figurename_multiwavelengths = "FloquetBandStructure-MultiWavelength-k" + str(
            kpoints) + "MPI-" + str(round(MPI_numbers[0], 0))+"-Field"+str(round(field,0))
            ax2.set_ylim(Materials[material]['ymin_ymax'])  # Si
            fig2.tight_layout()
            if(StorePlot):
                fig2.savefig(figurename_multiwavelengths + ".eps")
                fig2.savefig(figurename_multiwavelengths + ".png")
            if(ShowPlot):
                fig2.show()
            del fig2, ax2
        color_idx+=1

        

        # if(StorePlot):
        




        # print Header+"Preparing (omega,Efield) values that should disable the tunneling in the sample..."
        # plt.figure()

        print("=== Assessing tunneling stopping possibilities (Rabi oscillation) ... ===")
        Efield_AU_log = np.linspace(np.log10(Efield_AU_min), np.log10(Efield_AU_max), num_fields) #we need quasi-continous interval for RabiFrequencyMatrices
        Efield_AU = np.power(10., Efield_AU_log)

    if (PlotMatrixElements):
        figurename = str(material) + "-MatrixElements-k-" + str(kpoints)
        fig3, ax3 = plt.subplots()
        #ax3.matshow(matrixelements)
        im3 = ax3.imshow(np.log10(matrixelements))
        #plt.colorbar()
        fig3.colorbar(im3)
        plt.tight_layout()
        if(StorePlot):
            plt.savefig(figurename + ".eps")
            plt.savefig(figurename + ".png")

    if (PlotRabiFrequencyMatrices):
        for MPI_number in np.arange(-np.max(MPI_numbers), np.max(
                MPI_numbers) + 1):  # [Kato, PRB 100 081203(R) tells to use -MPI_numbers..MPI_numbers]
            tunneling_stops = []

            figurename = str(material) + "-MatrixElements-k-" + str(kpoints) + "MPI-" + str(round(MPI_number, 0))
            # for BesselRootsID in np.arange(1,BesselRootsID_max):
            BesselZeros = BesselJZeros(MPI_number, BesselRootsID_max)
            # print "BesselRoot #"+str(BesselRootsID_max)
            logger.debug("BesselJZeros(" + str(MPI_number) + ", " + str(BesselRootsID_max) + ")=" + str(BesselZeros))
            for RootNumber in np.arange(0, np.size(BesselZeros)):  # NOTE: we mix the various roots of Bessel function
                tunneling_stops.append(np.divide(matrixelements, BesselZeros[
                    RootNumber]))  # = omega**2/E. #NOTE: concatenation of all the roots of Bessel.
            logger.info(np.shape(tunneling_stops))
            # print tunneling_stops #omega**2/E (in at.u.)
            # plt.figure()
            plt.matshow(np.log10(tunneling_stops[MPI_number]))  # omega**2/E
            plt.title(r'$n=$' + str(MPI_number))  # +", slice "+str(BesselRootsID_max-1))
            plt.colorbar()
            plt.tight_layout()
            if(StorePlot):
                plt.savefig(figurename + ".eps")
                plt.savefig(figurename + ".png")

            figurename = str(material) + "-BesselTunnelingStop-" + str(kpoints) + "MPI-" + str(
                round(MPI_number, 0)) + "-BesselRoot-" + str(BesselRootsID_max)

            logger.info("Plotting couples of (E, omega) that blocks transition with " + str(MPI_number) + " photons.")
            # for wavelength_SI in wavelength_SI_t: #for each wavelength
            plt.figure()

            # Number of colors we need. 
            NumberOfElements = (matrix_side / 2)
            cm = get_cmap('hsv')
            colors = [cm(1. * i / (NumberOfElements)) for i in range(NumberOfElements)]
            xy = list(range(NumberOfElements))
            colorlist = [colors[x] for x in xy]

            plot1 = []
            for trans_i in np.arange(0, matrix_side / 2):  # transition from i
                for trans_j in np.arange(matrix_side / 2, matrix_side, 1):  # transition to j
                    for RootNumber in np.arange(0, np.size(BesselZeros)):
                        xx = Field_AU_to_SI(Efield_AU) * 1E-9
                        yy = au.Energy_Hartree_to_eV(np.sqrt(tunneling_stops[RootNumber][trans_i][trans_j] * Efield_AU))
                        plot_field = plt.scatter(xx, yy, s=SizeOfDots,
                                                 label=str(trans_i) + r"$\rightarrow$" + str(trans_j), c=colorlist[
                                trans_j - matrix_side / 2])  # colorlist[trans_j] is also possible #*matrix_side+trans_j]) #energies [eV]
                    plot1 += [plot_field]  # stores the last plot caption in a series of roots

            labels1 = [l.get_label() for l in plot1]
            plt.ylabel(r"$\hbar \omega$ (eV)")
            # plt.ylabel(r"$\lambda$ (m)")
            plt.legend(plot1, labels1, loc="best", ncol=2)
            plt.xlabel(r"$E$ (V/nm)")
            # plt.xscale('log')
            plt.ylim((0, 1.5))
            plt.xlim((0, 5))
            plt.title(r"Tunneling $(E,\omega)$, $n=$" + str(MPI_number))
            if(StorePlot):
                plt.savefig(figurename + ".eps")
                plt.savefig(figurename + ".png")
            # plt.show()


# }}}

## Plot the band structure at rest, along with the Floquet band structure for 1 value of field. 
## NOTE: Hypothesis that KpointPath.line_id = Kgrid.k_id. This can be wrong. Also check the K-point path in INPUT of octopus
# The best way to check is to plot bandstructure kx,ky,kz in space and compare with k-point path
# generated at https://tools.materialscloud.org/seekpath/compute/process_structure/
# @param field:         electric field amplitude (SI unit)
# @param wavelength_SI: laser wavelength (SI unit)
def MultipleKpointsOneField(field, wavelength_SI):  # {{{
    Header = "[plotFloquetBandStructure: MultipleKpointsOneField()] "
    global filename, listoffiles, Complex_valued_ME, unocc_states, Z_electrons, verbose, bandstructure, kpoints_manual, \
        KpointPath_len, Enable_A2, OpticalIndexFile
    CaptureDirectGap = False #We force it here since we need a number of k-point values for this function.

    # Capture all the k-point files for one direction
    try:
        kpoints_number_inp = libString.LookUpStringInFile("info", "Total number of k-points")
        kpoints_total_number = int(kpoints_number_inp.split("=")[1].strip())  # conversion to integer. Contains the total number of kpoints FOUND in the run.
        logger.info("** Info: total number of k-points = "+str(kpoints_total_number))
    except:
        logger.error("** Error: total number of k-points could not be detected.")
        exit()


    # Loading the bandstructure given by octopus (contains: kpointpath_index, kx, ky, kz, eigenvalues)
    bands = np.loadtxt(bandstructure, skiprows=1)
    logger.info((len(bands)))
    # bands_t=np.zeros()
    # exit()
    # 1. Build list of kpoint-index employed by octopus in <info> from the <bandstructure>.

    # It seems kpointpath is actually available in info: it is just located at the end of the last kpoints_number_array. 

    # kpointpath_index=bands[0,:]
    # kpointpath_kx   =bands[1,:]
    # kpointpath_ky   =bands[2,:]
    # kpointpath_kz   =bands[3,:]

    # bands_t = np.array(kpointpath_index, kpointpath_kid, kpointpath_kx, kpointpath_ky, kpointpath_kz, bands[4:,:]) #just adding a column to the object
    # 2. Plot the bandstructure

    # print("kpoints_number"+str(kpoints_number))^
    # kpoints_number_array = kpoints_number_inp #.split("_")
    # print(kpoints_number_array)
    # exit()
    # kpoints_number_array = kpoints_number_array[0].split("k")
    # kpoints       = int(kpoints_number_array[0])
    # print(kpoints)
    # except:
    # kpoints                = kpoints_manual #manual input of k-point grid unit (assumes same in 3 directions)
    # print Header+"** Error: failed to capture the number of k points automatically."
    # exit()

    # KpointPath_list: we know it is located after the 3D kgrid made kpoints**3 entries (e.g. 16**3 = 4096 points).

    # Compute the number of 3D k-point files.
    try:
        kpoints3D_number = libString.CaptureKPointGridSizeFromInfoFile("info")
    except:
         try:
             kpoints3D_number=libString.CaptureKPointGridSizeFromINP("../inp")
         except: #if file was not found
             logger.warning("** Warning: could not capture size of K-point grid in INP file. Most probably absent.")
             kpoints3D_number = kpoints_manual
             logger.warning("**          Employing manual k-point number as last resort: ")
    logger.info("** Number of 3D kpoints: "+str(kpoints3D_number))

    # Capture the total number of k-points (available in info file)
    ## Captures the number of 3D k-points in the 2D-path
    # Hypothesis: the path starts with Gamma point. Lets capture the kpt_id when finding 0|0|0 in the list of info

    # Deduce the number of points in the k-point path.
    KpointPath_len = kpoints_total_number - kpoints3D_number**3
    KpointPath_list = np.arange(kpoints3D_number ** 3 + 1,
                                kpoints_total_number) #NOTE: this assumes no other path is added at the end!
                                #kpoints3D_number ** 3 + KpointPath_len) #TODO: KpointPath_len should be fetched automatically!!
    # NOTE: this may wrongly assume a matching between Kgrid.k_id and KpointPath.line_id. I checked position of repeated Gamma points in the path. Dipolar elements were stricly equal at these positions. 
    # breakpoint()
    logger.info(("KPointPath_list: " + str(len(KpointPath_list))))
    ## Obtain a list of all dipolar data for a single direction for a list of k-points
    logger.info("Info: Capturing the files for the transition at given k-points.")
    k_point_files, k_point_files_idx = libString.ExtractFileForGivenKpoints(filename, kpoints3D_number, listoffiles,
                                                                        KpointPath_list, CaptureDirectGap, verbose,
                                                                        field_polarization_dir)

    # k_point_files contains a table of names for all k-point files in a single direction.
    # k_point_files_idx contains a table of indices matching with the k_points_files. Only selection in X direction was made.

    # We can also prepare a figure of the k-point position in 3D. Then we'll be able to check if it corresponds to
    # Fermi sphere from books.

    plt.figure()
    ax1 = plt.subplot(211)
    ax1.set_title("Band structure")
    ax1.set_ylabel("Energy (eV)")
    
    ax2 = plt.subplot(212)
    ax2.set_title("Non-equilibrium band structure")
    ax2.set_ylabel("Energy (eV)")
    ax2.set_xlabel("k-point path")

    # ax2.set_xlabel("k-point id")

    ## NOTE: if we assume that laser light only REMOVES some high-symmetry points from the ones of the space group, 
    ## NOTE: then we have to get the k-points from the file <bandstructure> since interpolation is needed between the points of high-symmetry). INP file does not provide this. 

    ## Plotting the band structure
    # print bands[0,4:] #shows all eigenvalues of 1st line
    # print bands[:,0] #shows all kpoints norm
    # exit()

    ## THIS IS the proper BAND STRUCTURE from <bandstructure>. We just plot it. 
    #plt.title("Band structure")
    # ax1.set_xticklabels([])
    #ax1.set_ylabel("Energy (eV)")
    for eig in np.arange(4, len(bands[0, :])):
        ax1.scatter(bands[:, 0], au.Energy_Hartree_to_eV(bands[:, eig]),
                    s=2)  # each band will be plotted with a different color
    # plt.show()

    ## TO ENSURE THAT KPointPath is well located at the end of k_id in <info>, we plot them both. 

    # Build a table that associates the path of k-point (found in inp or bandstructure file) with index of the corresponding k-points.
    ## Index matching will be done by comparing reduced kx,ky,kz between info and bandstructure. 

    # for each kpointpath k-point, we need to search the index that matches with <info_kx, info_ky, info_kz>
    # for kpt_id in np.arange(0,kpointpath_kx):
    # if(kpointpath_kx==info_kx and kpointpath_ky==info_ky and kpointpath_kz==info_kz):

    ### We want to keep the k-point path defined in bandstructure. 
    ### We parse bandstructure. 
    ### For each kpoint_reduced coords, we look for these coordinates in info, and associate it with k_point_index. 
    ### Then we simply change the loop below to work on this list of k-point indices. We will construct the Floquet bandstructure quite easily then. 

    ## NOTE: In case laser would change the high-symmetry points, we would need to visualize all the k-points. We need to access dipolar matrices anyway, not given by the bandstructure file. 
    idx_t = []
    # print(np.shape(idx_t))
    # exit()
    logger.debug(("** Debug: size of k_point_files: " + str(len(k_point_files))))
    logger.debug(("         size of k_point_files_idx: " + str(len(k_point_files_idx))))
    # print("         content of k_point_files: "+str(k_point_files))
    # print("         content of k_point_files_idx: "+str(k_point_files_idx))

    ## DOUBLE PARSING: 
    ## - k_idx belongs to the manifold of KPathPoints (sync with bandstructure)
    ## - k_point_files[k_idx] belongs to KPointGrid (sync with info and ks_me_dipole*)
    for k_idx in np.arange(0, len(k_point_files)):  # fake index for parsing the KPointPath
        # for k_point_file in k_point_files: #going through each kpoint. We here restrict the Floquet computation to the path of k-points.
        if(verbose==1):
            logger.debug(("** Info: k_idx in path: " + str(k_idx)))
            logger.debug(("** Info: k_point_file: " + str(k_point_files[k_idx])))
            logger.debug(("** Info: fetching matrix elements from " + str(k_point_files[k_idx])))
            logger.debug(("** Info: k_idx in full list is: " + str(k_point_files_idx[k_idx])))

        # print(Z_electrons, unocc_states, k_point_file, Complex_valued_ME)

        ## 5. Extract eigen values at k-point <kpoints_index> from "info" file
        if(verbose==1):
            logger.debug(("** Fetching eigenvalues for k-point " + str(k_point_files[k_idx])))
        eigenvalues = libString.GetEigenValue(k_point_files_idx[k_idx], Z_electrons + unocc_states, "info", verbose)
        if(verbose==1):
            logger.debug(eigenvalues)

        # Plot the band structure at rest using provided k-point path
        if(verbose==1):
            logger.debug(("** Info: shapes: " + str(np.shape(idx_t)) + ", " + str(np.shape(eigenvalues))))
        # idx_t=eigenvalues
        if(verbose==1):
            logger.debug((str(np.shape(idx_t)) + ", " + str(np.shape(eigenvalues))))

        # X: we shall plot as function of bands[0,:]. 
        #ax2.scatter(bands[k_idx, 0] * np.ones(len(eigenvalues)), au.Energy_Hartree_to_eV(eigenvalues), c="blue", s=0.1, label="Eq.")
        if (not PlotFloquetBandStructure):
            continue

        ## Loading matrix elements for the k-point.
        if(verbose==1):
            logger.info(("** Info: fetching matrix elements."))
        matrixelements, matrix_side = libString.ReadDipolarMatrixElements(Z_electrons, unocc_states, [k_point_files[k_idx]],
                                                                Complex_valued_ME)
        # Compute the Floquet band structure for each selected k-point file, but only for one field
        # wavelength_SI = SelectedWavelengths
        epsilon_Si0 = libOpticalData.ImportPalikDatabase_epsilon_fromNK(SelectedWavelengths,
                                                                        OpticalIndexFile)
        ## WARNING: replacing by LDA is DIFFICULT, since PRL Rubio is defined down to 2 eV. Nothing below!! #TODO: this should be material-dependent. Indicate path to fetch the right file in library at the top.
        RefractiveIndex_t = np.real(np.sqrt(epsilon_Si0))  # This is a continuum formulation from Palik data
        RefractiveIndex = RefractiveIndex_t #OpticalIndex[wavelength_SI]
        if (UseRefractiveIndex):
            Efield_SI = field / np.real(np.sqrt(epsilon_Si0)) # np.sqrt(OpticalIndex[wavelength_SI])
        else:
            Efield_SI = field
        # for MPI_number in MPI_numbers:
        MPI_number = MPI_numbers[0]
        eigenvalues, replicas, dipoles = ComputeFloquetBandStructure(filename, nb_atoms, Z_electrons,
                                                                     k_point_files_idx[k_idx], unocc_states,
                                                                     MPI_number, Efield_SI, wavelength_SI,
                                                                     num_time_steps, eigenvalues, matrixelements,
                                                                     Complex_valued_ME, field_polarization_dir,
                                                                     spin_occupation, verbose, Enable_A2)

        if (PlotIntensity):
            unitx = 1E-4
            Efield_SI = 0.5 * c * epsilon_0 * Efield_SI ** 2 * RefractiveIndex
        else:
            unitx = 1E0
            Efield_SI = Efield_SI * np.sqrt(RefractiveIndex)
        ax2.scatter(np.ones(np.size(replicas)) * bands[k_idx, 0], au.Energy_Hartree_to_eV(np.real(replicas)), c="red",
                    s=0.1*SizeOfDots, label="Floquet replicates")  # these are the interesting lines for the final plot
        #breakpoint()
    ax1.set_ylim(Materials[material]['ymin_ymax'])
    ax2.set_ylim(Materials[material]['ymin_ymax'])
    ax1.set_xticks(np.linspace(0., 1., len(K_labels)))
    ax2.set_xticks(np.linspace(0., 1., len(K_labels)))
    ax1.set_xticklabels(K_labels, fontsize='15')
    ax2.set_xticklabels(K_labels, fontsize='15')
    plt.tight_layout()
    figurename = "FloquetBandStructureComplete-" + str(round(wavelength_SI * 1E9, 0)) + "nm-k" + str(
        kpoints3D_number) + "MPI-" + str(round(MPI_number, 0)) + "-Field" + str(1E-9 * round(field, 0)) + "VperNm"
    if (not LogPlotting):
        figurename += "NoLog"
    if (not PlotIntensity):
        figurename += "-Field"
    #plt.savefig(figurename + ".eps")
    if(StorePlot):
        plt.savefig(figurename + ".png")
    if(ShowPlot):
        plt.show()


# }}}

if(__name__=="__main__"):  #runs only if file is called directly
    if (SingleKpoint == 1):
        SingleKpointManyFields()
    else:
        for wavelength in [1030e-9]: #SelectedWavelengths:
            # MultipleKpointsOneField(1E8)
            MultipleKpointsOneField(1E9, wavelength)
            # MultipleKpointsOneField(5E9)
            MultipleKpointsOneField(1E10, wavelength)
            MultipleKpointsOneField(1E11, wavelength)

