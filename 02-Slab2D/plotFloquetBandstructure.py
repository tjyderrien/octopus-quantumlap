#!/usr/bin/env python3
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
# Provides routines and prepare dressed band structure from ground state matrix elements of a solid crystal

from matplotlib.cm import get_cmap
from scipy.constants import c, epsilon_0, h
from scipy.special import jn_zeros as BesselJZeros

from importlib import import_module
IntensityToField = import_module("spp-extended-theory.Keldysh.libKeldyshZhukov.IntensityToField")
from libPlotting import *
from libPlottingSettings import *
ComputeFloquetBandStructure = import_module('spp-extended-theory.Keldysh.libStark.ComputeFloquetBandStructure')
from libString import *

PlotMatrixElements = True
PlotFloquetBandStructure = True  # Provide field dependent band structure at Gamma Point
PlotRabiFrequencyMatrices = False  # for studying Wannier-Stark localization, i.e. suppression of tunneling with zeros of the Bessel function

# Filenames = "ks_me_multipoles.*"
Filenames = "ks_me_dipole.*"
bandstructure = "bandstructure"  # NOTE: file may be absent sometimes
LocalFolder = "./"
filename = Filenames.split('.')[0]
# print Filenames
# 1. Capture the static/dipole files. 
listoffiles = GetOrderedListOfFiles(Filenames)
del Filenames
# print listoffiles

"""
Explanation of the number of method employed by Octopus for indexing the k-grid files:
- Numbering of the files is N_files = (kgrid) + (kpoint path)
- 1,2,3 directions for dipolar matrix elements (x,y,z), added a suffix of the dipolar files ".x###". 
- Example: k-grid size=(16 x 16 x 16) = 4096 files for 3D k-band structure. 
- One k-point path is given at the end of the numbering. e.g. 1001 files remain for k-point path. 
"""

## Building materials dictionnary
dirname=os.path.realpath(".")

print(("** Info: directory name: "+str(dirname)))

if "03-Bulk3D" in dirname:
    material = "Si"
elif "09-Au" in dirname: 
    material = "Au"
elif "06-Mo" in dirname: 
    material = "Mo"
elif "14-SiO2" in dirname: 
    material = "SiO2"
else: 
    print(Header+"Which material are you trying to analyze?")
    
# material = "Au"  # SiO2 | Si

# exit()
wavelength_SI_t = [800e-9]  # , 400e-9] #484e-9, 343e-9] #10e-6

num_time_steps = 100 #should be small enough to ensure convergence
num_fields = 300  # number of electric field values for the plotting

Efield_SI_min = IntensityToField(1E10 * 1E4)
Efield_SI_max = IntensityToField(1E19 * 1E4)
Efield_AU_min = Field_SI_to_AU(Efield_SI_min)
Efield_AU_max = Field_SI_to_AU(Efield_SI_max)

Complex_valued_ME = False #True for non-Hermitian eigenvalues?
field_polarization_dir = 0  # 0: x, 1: y, 2: z
spin_occupation = 2  # 1: levels with same spin should split, 2: levels of different spins don't split.

Renormalize = False  # True: limits the plots maximum as function of wavelength
PlotIntensity = True  # True: plot with I. False: plot with E.
SizeOfDots = 0.5  # purely cosmetic!
LogPlotting = 1
verbose = 0  # Show info and debug messages.
SingleKpoint = 1  # 0: plotting the whole band structure for a certain field. 1: plot one single kpoint for many fields.

BesselRootsID_max = 10

########### DONT TOUCH BELOW, thanks! #############

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
Si = {'MPI_numbers': [1, 2, 3], 'kpoints_manual': 16,
      'nb_atoms': 1, 'unocc_states': 4,
      'Z_electrons': 4, 'OpticalIndex': OpticalIndex_Si,
      'ymin_ymax': (3, 8), 'CaptureDirectGap': True, 'KpointPath_len': 301}

SiO2 = {'MPI_numbers': [1, 2], 'kpoints_manual': 16,
        'nb_atoms': 1, 'unocc_states': 6,
        'Z_electrons': 24, 'OpticalIndex': OpticalIndex_SiO2,
        'ymin_ymax': (-8, 16), 'CaptureDirectGap': True, 'KpointPath_len': 301}

# In metals, we are interested in plotting band structure on the full path of high symmetry k-points
Mo = {'MPI_numbers': [1], 'kpoints_manual': 16,
      'nb_atoms': 1, 'unocc_states': 3,
      'Z_electrons': 8, 'OpticalIndex': OpticalIndex_Mo,
      'ymin_ymax': (-30, 30), 'CaptureDirectGap': False, 'KpointPath_len': 601,
      'K_labels': [r'$\Gamma$', r'$H$', r'$N$', r'$\Gamma$', r'$P$', r'$H|P$', r'$N$']}

Au = {'MPI_numbers': [1], 'kpoints_manual': 16,
      'nb_atoms': 1, 'unocc_states': 3,
      'Z_electrons': 14 - 3, 'OpticalIndex': OpticalIndex_Au,
      'ymin_ymax': (-2, 35), 'CaptureDirectGap': False, 'KpointPath_len': 1001,
      'K_labels': [r'$\Gamma$', r'$X$', r'$W$', r'$U|K$', r'$\Gamma$', r'$L$', r'$U$', r'$W$', r'$L$', r'$K|U$',
                   r'$X$']}

Materials = {'Si': Si, 'SiO2': SiO2, 'Mo': Mo, 'Au': Au}  # this contains all materials properties we need to automatize

# 2. Check consistency of the number of files with k-point employed in the computation
nb_atoms = Materials[material]['nb_atoms']  # 1 #number of atoms #we actually treat a bulk, so it's one.
Z_electrons = Materials[material][
    'Z_electrons']  # SiO2: 24 electrons. #Si: 4 #number of modeled valence electrons in the simulation box
# (3p)² x 2 for Si (biatomic) lattice. 
# SiO2 ? matrix is 30 e- states, 6 unocc. 
unocc_states = Materials[material][
    'unocc_states']  # 4 #SiO2: 6. #Si: 4. This is defined in the inp. If not, then it's zero.

CaptureDirectGap = Materials[material]['CaptureDirectGap']

MPI_numbers = Materials[material]['MPI_numbers']  # [1, 2, 3]
kpoints_manual = Materials[material]['kpoints_manual']  # 16 #number of kpoints in the grid
OpticalIndex = Materials[material]['OpticalIndex']
KpointPath_len = Materials[material]['KpointPath_len']
# print
K_labels = Materials[material]['K_labels']


## Prepare energy diagram as function of laser field strength at a single K point
def SingleKpointManyFields():  # {{{
    global Efield_SI_min, Efield_SI_max, kpoints_manual, Z_electrons, unocc_states
    ## Capturing the Gamma k-point file in results of Octopus 
    # TODO: generalize to a desired k-point
    try:
        kpoints_number_inp = LookUpStringInFile("info",
                                                "Dimensions of the k-point grid")  # WARNING: it seems this has changed
        for i in np.arange(0, 10):  # TODO: replace with .strip() ?
            kpoints_number_inp = kpoints_number_inp.replace("  ", " ")
            if (verbose == 1):
                print(Header + "** Info: decreasing number of spaces")
                print(kpoints_number_inp)
        kpoints_number_array = kpoints_number_inp.split("=")
        kpoints_number_array = kpoints_number_array[1].split(" ")
        kpoints = int(kpoints_number_array[1])
        print(kpoints)
    except:
        kpoints = kpoints_manual  # manual input of number of k-grid points in 1 direction (assumes same resolution in 3 directions)
        print(
            Header + "** Error: failed to capture the number of k points automatically. Manual input is employed: kpoints=" + str(
                kpoints))
        # exit()

    print("Info: Capturing the file for the transition of interest.")
    gamma_point_files, gamma_point_files_index = ExtractFileForGivenKpoints(filename, kpoints, listoffiles, [1],
                                                                            CaptureDirectGap, verbose,
                                                                            field_polarization_dir)

    matrixelements, matrix_side = ReadDipolarMatrixElements(Z_electrons, unocc_states, gamma_point_files,
                                                            Complex_valued_ME)

    ## 5. Extract eigen values at k-point <kpoints_index> from "info" file
    eigenvalues = GetEigenValue(gamma_point_files_index[0], Z_electrons + unocc_states, "info", verbose)
    if (verbose == 1):
        print(eigenvalues)

    # tunneling_stops=[]
    if (PlotFloquetBandStructure):
        for wavelength_SI in wavelength_SI_t:

            RefractiveIndex = OpticalIndex[wavelength_SI]

            ## Renormalize the intensity range using refractive index
            if (Renormalize):
                Efield_SI_max = Efield_SI_max / np.sqrt(OpticalIndex[wavelength_SI])

            if (LogPlotting):
                Efield_SI_log = np.linspace(np.log10(Efield_SI_min), np.log10(Efield_SI_max), num_fields)
                Efield_SI = np.power(10., Efield_SI_log)
            else:
                Efield_SI = np.linspace(Efield_SI_min, Efield_SI_max, num_fields)

            ## Plot the splitted energy levels at Gamma.
            for MPI_number in MPI_numbers:
                plt.figure()
                for field in Efield_SI:
                    eigenvalues, replicas, dipoles = ComputeFloquetBandStructure(filename, nb_atoms, Z_electrons,
                                                                                 kpoints, unocc_states, MPI_number,
                                                                                 field, wavelength_SI, num_time_steps,
                                                                                 eigenvalues, matrixelements,
                                                                                 Complex_valued_ME,
                                                                                 field_polarization_dir,
                                                                                 spin_occupation, verbose)
                    # Plotting using external intensity, whence internal intensity acts on the sample
                    if (PlotIntensity):
                        unitx = 1E-4
                        field = 0.5 * c * epsilon_0 * field ** 2 * RefractiveIndex
                    else:
                        unitx = 1E0
                        field = field * np.sqrt(RefractiveIndex)
                    plt.scatter(np.ones(np.size(replicas)) * field * unitx, Energy_Hartree_to_eV(np.real(replicas)),
                                c="red", s=SizeOfDots,
                                label="Floquet replicates")  # these are the interesting lines for the final plot

                    # TODO: Prepare the table [wavelength, n, E, Rabi Freq such as BesselJ(n,RabiFreqs_SI/omegaLaser)=0. #BruteForce method. 
                    # TODO: Do we need the function BesselJZeros(n,x)? To search omegaLaser that induce a stop in tunneling! 
                    # tunneling_stops.append(np.sqrt(dipoles * field_AU / BesselJZeros(n,k)))

                if (PlotIntensity):
                    Intensity_min = 0.5 * c * epsilon_0 * Efield_SI[0] ** 2 * RefractiveIndex
                else:
                    Intensity_min = Efield_SI[0] * np.sqrt(RefractiveIndex)

                # Plotting eigenvalues at rest
                plt.scatter(np.ones(np.size(eigenvalues)) * Intensity_min * unitx,
                            Energy_Hartree_to_eV(np.real(eigenvalues)), c="black", s=10, label="GS eigenvalues")

                if (LogPlotting):
                    plt.xscale('log')

                plt.title(r"$n=$" + str(round(MPI_number, 0)) + r", $\lambda=$" + str(
                    int(wavelength_SI * 1E9)) + " nm, k=" + str(kpoints) + r"$^3$")
                figurename = "FloquetBandStructure-" + str(round(wavelength_SI * 1E9, 0)) + "nm-k" + str(
                    kpoints) + "MPI-" + str(round(MPI_number, 0))
                if (not LogPlotting):
                    figurename += "NoLog"
                if (not PlotIntensity):
                    figurename += "-Field"
                    plt.xlabel("Field (V/m)")
                else:
                    plt.xlabel(r"Laser intensity (W/cm$^{2}$)")
                plt.ylabel(r"Energy levels (eV)")
                plt.ylim(Materials[material]['ymin_ymax'])  # Si
                # plt.ylim((-8,8)) #SiO2
                plt.tight_layout()
                plt.savefig(figurename + ".eps")
                plt.savefig(figurename + ".png")
                # Plotting dressed levels directly on the figure, from n=1 to MPI_number
                for MPI in np.arange(1, MPI_number + 1):
                    plt.scatter(np.ones(np.size(eigenvalues)) * Intensity_min * unitx,
                                Energy_Hartree_to_eV(np.real(eigenvalues)) + MPI * h * c / wavelength_SI / e, c="red",
                                s=10, label="non-interacting replicates")
                    plt.scatter(np.ones(np.size(eigenvalues)) * Intensity_min * unitx,
                                Energy_Hartree_to_eV(np.real(eigenvalues)) - MPI * h * c / wavelength_SI / e, c="blue",
                                s=10, label="non-interacting replicates")

                plt.savefig(figurename + "-Suppl" + ".eps")
                plt.savefig(figurename + "-Suppl" + ".png")
                # plt.legend()
            # plt.show()
        # print Header+"Preparing (omega,Efield) values that should disable the tunneling in the sample..."
        # plt.figure()

        print("=== Assessing tunneling stopping possibilities... ===")
        Efield_AU_log = np.linspace(np.log10(Efield_AU_min), np.log10(Efield_AU_max), num_fields)
        Efield_AU = np.power(10., Efield_AU_log)

    if (PlotMatrixElements):
        figurename = str(material) + "-MatrixElements-k-" + str(kpoints)
        plt.figure()
        plt.matshow(matrixelements)
        plt.colorbar()
        plt.tight_layout()
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
            print("BesselJZeros(" + str(MPI_number) + ", " + str(BesselRootsID_max) + ")=" + str(BesselZeros))
            for RootNumber in np.arange(0, np.size(BesselZeros)):  # NOTE: we mix the various roots of Bessel function
                tunneling_stops.append(np.divide(matrixelements, BesselZeros[
                    RootNumber]))  # = omega**2/E. #NOTE: concatenation of all the roots of Bessel.
            print(np.shape(tunneling_stops))
            # print tunneling_stops #omega**2/E (in at.u.)
            # plt.figure()
            plt.matshow(np.log10(tunneling_stops[MPI_number]))  # omega**2/E
            plt.title(r'$n=$' + str(MPI_number))  # +", slice "+str(BesselRootsID_max-1))
            plt.colorbar()
            plt.tight_layout()
            plt.savefig(figurename + ".eps")
            plt.savefig(figurename + ".png")

            figurename = str(material) + "-BesselTunnelingStop-" + str(kpoints) + "MPI-" + str(
                round(MPI_number, 0)) + "-BesselRoot-" + str(BesselRootsID_max)

            print("Plotting couples of (E, omega) that blocks transition with " + str(MPI_number) + " photons.")
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
                        yy = Energy_Hartree_to_eV(np.sqrt(tunneling_stops[RootNumber][trans_i][trans_j] * Efield_AU))
                        plot_field = plt.scatter(xx, yy, s=10 * SizeOfDots,
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
            plt.savefig(figurename + ".eps")
            plt.savefig(figurename + ".png")
            # plt.show()


# }}}

## Plot the band structure at rest, along with the Floquet band structure for 1 value of field. 
## BUG: Hypothesis that KpointPath.line_id = Kgrid.k_id is wrong. TODO: Find where this is employed in the algorithm. 
def MultipleKpointsOneField(field):  # {{{
    Header = "[plotFloquetBandStructure: MultipleKpointsOneField()] "
    global filename, listoffiles, Complex_valued_ME, unocc_states, Z_electrons, verbose, bandstructure, kpoints_manual, KpointPath_len
    # Capture all the k-point files for one direction
    # try:
    kpoints_number_inp = LookUpStringInFile("info", "Total number of k-points")
    kpoints_number = int(kpoints_number_inp.split("=")[1].strip())  # integer

    # Loading the bandstructure given by octopus (contains: kpointpath_index, kx, ky, kz, eigenvalues)
    bands = np.loadtxt(bandstructure, skiprows=1)
    print((len(bands)))
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

    # KpointPath_list: we know it is located after the kgrid description. No need for searching it. More efficient than searching algorithm with N**2 complexity. 
    KpointPath_list = np.arange(kpoints_manual ** 3 + 1,
                                kpoints_manual ** 3 + KpointPath_len + 1)  # WARNING: this may wrongly assume a matching between Kgrid.k_id and KpointPath.line_id
    print(("KPointPath_list: " + str(len(KpointPath_list))))
    ## Obtain a list of all dipolar data for a single direction for a list of k-points
    print("Info: Capturing the files for the transition at given k-points.")
    k_point_files, gamma_point_files_index = ExtractFileForGivenKpoints(filename, kpoints_manual, listoffiles,
                                                                        KpointPath_list, CaptureDirectGap, verbose,
                                                                        field_polarization_dir)

    # k_point_files contains a table of names for all k-point files in a single direction.
    # gamma_point_files_index contains a table of indices matching with the k_points_files. Only selection in X direction was made.

    plt.figure()
    ax1 = plt.subplot(212)
    plt.title("Full unsorted band structure")

    # ax1.set_xlabel("k-point id")
    ax1.set_ylabel("Energy (eV)")

    ax2 = plt.subplot(211)
    plt.title("Band structure")

    ## NOTE: if we assume that laser light only REMOVES some high-symmetry points from the ones of the space group, 
    ## NOTE: then we have to get the k-points from the file <bandstructure> since interpolation is needed between the points of high-symmetry). INP file does not provide this. 

    ## Plotting the band structure
    # print bands[0,4:] #shows all eigenvalues of 1st line
    # print bands[:,0] #shows all kpoints norm
    # exit()

    ## THIS IS the proper BAND STRUCTURE from <bandstructure>. We just plot it. 
    plt.title("Band structure")
    # ax2.set_xticklabels([])
    ax2.set_xlabel("k-point path")
    ax2.set_ylabel("Energy (eV)")
    for eig in np.arange(4, len(bands[0, :])):
        ax2.scatter(bands[:, 0], Energy_Hartree_to_eV(bands[:, eig]),
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
    print((Header + "** Info: size of k_point_files: " + str(len(k_point_files))))
    print((Header + "         size of gamma_point_files_index: " + str(len(gamma_point_files_index))))
    # print(Header+"         content of k_point_files: "+str(k_point_files))
    # print(Header+"         content of gamma_point_files_index: "+str(gamma_point_files_index))

    ## DOUBLE PARSING: 
    ## - k_idx belongs to the manifold of KPathPoints (sync with bandstructure)
    ## - k_point_files[k_idx] belongs to KPointGrid (sync with info and ks_me_dipole*)
    for k_idx in np.arange(0, len(k_point_files)):  # fake index for parsing the KPointPath
        # for k_point_file in k_point_files: #going through each kpoint. We here restrict the Floquet computation to the path of k-points.
        print((Header + "** Info: k_idx in path: " + str(k_idx)))
        print((Header + "** Info: k_point_file: " + str(k_point_files[k_idx])))
        print(("** Info: fetching matrix elements from " + str(k_point_files[k_idx])))
        print(("** Info: k_idx in full list is: " + str(gamma_point_files_index[k_idx])))

        # print(Z_electrons, unocc_states, k_point_file, Complex_valued_ME)

        ## 5. Extract eigen values at k-point <kpoints_index> from "info" file
        print(("** Fetching eigenvalues for k-point " + str(k_point_files[k_idx])))
        eigenvalues = GetEigenValue(gamma_point_files_index[k_idx], Z_electrons + unocc_states, "info", verbose)
        # if(verbose==1):
        print(eigenvalues)

        # Plot the band structure at rest using provided k-point path
        print((Header + "** Info: shapes: " + str(np.shape(idx_t)) + ", " + str(np.shape(eigenvalues))))
        # idx_t=eigenvalues

        print((Header + str(np.shape(idx_t)) + ", " + str(np.shape(eigenvalues))))

        # X: we shall plot as function of bands[0,:]. 
        ax1.scatter(bands[k_idx, 0] * np.ones(len(eigenvalues)), Energy_Hartree_to_eV(eigenvalues), c="blue", s=0.1)
        if (not PlotFloquetBandStructure):
            continue

        ## Loading matrix elements for the k-point.
        print((Header + "** Info: fetching matrix elements."))
        matrixelements, matrix_side = ReadDipolarMatrixElements(Z_electrons, unocc_states, [k_point_files[k_idx]],
                                                                Complex_valued_ME)
        # Compute the Floquet band structure for each selected k-point file, but only for one field
        wavelength_SI = wavelength_SI_t[0]
        RefractiveIndex = OpticalIndex[wavelength_SI]
        if (Renormalize):
            Efield_SI = field / np.sqrt(OpticalIndex[wavelength_SI])
        else:
            Efield_SI = field
        # for MPI_number in MPI_numbers:
        MPI_number = MPI_numbers[0]
        eigenvalues, replicas, dipoles = ComputeFloquetBandStructure(filename, nb_atoms, Z_electrons,
                                                                     gamma_point_files_index[k_idx], unocc_states,
                                                                     MPI_number, Efield_SI, wavelength_SI,
                                                                     num_time_steps, eigenvalues, matrixelements,
                                                                     Complex_valued_ME, field_polarization_dir,
                                                                     spin_occupation, verbose)

        if (PlotIntensity):
            unitx = 1E-4
            Efield_SI = 0.5 * c * epsilon_0 * Efield_SI ** 2 * RefractiveIndex
        else:
            unitx = 1E0
            Efield_SI = Efield_SI * np.sqrt(RefractiveIndex)
        ax1.scatter(np.ones(np.size(replicas)) * bands[k_idx, 0], Energy_Hartree_to_eV(np.real(replicas)), c="red",
                    s=SizeOfDots, label="Floquet replicates")  # these are the interesting lines for the final plot
    ax1.set_ylim(Materials[material]['ymin_ymax'])
    ax2.set_ylim(Materials[material]['ymin_ymax'])
    ax1.set_xticks(np.linspace(0., 1., len(K_labels)))
    ax2.set_xticks(np.linspace(0., 1., len(K_labels)))
    ax1.set_xticklabels(K_labels, fontsize='15')
    ax2.set_xticklabels(K_labels, fontsize='15')
    plt.tight_layout()
    figurename = "FloquetBandStructureComplete-" + str(round(wavelength_SI * 1E9, 0)) + "nm-k" + str(
        kpoints_manual) + "MPI-" + str(round(MPI_number, 0)) + "-Field" + str(1E-9 * round(field, 0)) + "VperNm"
    if (not LogPlotting):
        figurename += "NoLog"
    if (not PlotIntensity):
        figurename += "-Field"
    plt.savefig(figurename + ".eps")
    plt.savefig(figurename + ".png")
    plt.show()


# }}}

if (SingleKpoint == 1):
    SingleKpointManyFields()
else:
    # MultipleKpointsOneField(1E8)
    MultipleKpointsOneField(1E9)
    # MultipleKpointsOneField(5E9)
    MultipleKpointsOneField(1E10)
    MultipleKpointsOneField(1E11)

