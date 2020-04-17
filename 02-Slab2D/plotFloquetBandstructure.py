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

## @package plotFloquetBandstructure
# Provides routines and prepare dressed band structure from ground state matrix elements

import numpy as np
import glob
import numpy.linalg as LA

from libAtomicUnits      import *
from libPlotting         import *
from libPlottingSettings import *
from libString           import *
from libKeldyshZhukov    import IntensityToField
from libStark            import Kronecker, ComputeFloquetBandStructure

#Filenames = "ks_me_multipoles.*"
Filenames = "ks_me_dipole.*"
LocalFolder="./"
filename      = Filenames.split('.')[0]
print Filenames
# 1. Capture the static/dipole files. 
listoffiles=GetOrderedListOfFiles(Filenames)
print listoffiles

## Building materials dictionnary
material="SiO2" #SiO2 | Si

#exit()
wavelength_SI_t        = [3200e-9, 1600e-9, 800e-9, 484e-9, 400e-9, 343e-9]
#wavelength_SI_t        = [484e-9] #, 400e-9, 343e-9]

num_time_steps         = 100
num_fields             = 300 #number of electric field values for the plotting

Efield_SI_min          = IntensityToField(1E10*1E4)
Efield_SI_max          = IntensityToField(1E15*1E4)

Complex_valued_ME      = False
field_polarization_dir = 0 #0: x, 1: y, 2: z
spin_occupation        = 2 #1: levels with same spin should split, 2: levels of different spins don't split.

SizeOfDots             = 0.5   #purely cosmetic!
LogPlotting            = 1
verbose                = 0 #Show info and debug messages. 

########### DONT TOUCH BELOW, thanks! #############

OpticalIndex_Si={ #NOTE: these are experimental data for Si (from Palik). They should be from linear TDDFT response !
    343e-9: 5.2749, 400e-9: 5.5674, 483e-9: 4.4028, 484e-9: 4.3964,
    800e-9: 3.6924167231905, 1600e-9: 3.4826, 
    2200e-9: 3.4548, 2600e-9: 3.4457, 3200e-9: 3.43438589865999}

OpticalIndex_SiO2={ #NOTE: these are experimental data for SiO2 (from refractiveindex.info). They should be from linear TDDFT response !
    343e-9: 1.5669, 400e-9: 1.5577, 483e-9: 1.5499, 484e-9: 1.5498,
    800e-9: 1.5383, 1600e-9: 1.5270, 
    2200e-9: 1.5175, 2600e-9: 1.5094, 3200e-9: 1.4996}

Si={'MPI_numbers': [1, 2, 3], 'kpoints_manual': 16, 
    'nb_atoms': 1, 'unocc_states': 4, 
    'Z_electrons': 4, 'OpticalIndex': OpticalIndex_Si, 
    'ymin_ymax': (3,8)}

SiO2={'MPI_numbers': [1], 'kpoints_manual': 8, 
    'nb_atoms': 1, 'unocc_states': 6, 
    'Z_electrons': 24, 'OpticalIndex': OpticalIndex_SiO2, 
    'ymin_ymax': (-8,8)}

Materials={'Si': Si, 'SiO2': SiO2} #this contains all materials properties we need to automatize

# 2. Check consistency of the number of files with k-point employed in the computation
nb_atoms      = Materials[material]['nb_atoms'] #1 #number of atoms #we actually treat a bulk, so it's one.
Z_electrons   = Materials[material]['Z_electrons'] #SiO2: 24 electrons. #Si: 4 #number of modeled valence electrons in the simulation box 
# (3p)² x 2 for Si (biatomic) lattice. 
# SiO2 ? matrix is 30 e- states, 6 unocc. 
unocc_states           = Materials[material]['unocc_states'] #4 #SiO2: 6. #Si: 4. This is defined in the inp. If not, then it's zero. 

MPI_numbers             = Materials[material]['MPI_numbers'] #[1, 2, 3]
kpoints_manual          = Materials[material]['kpoints_manual'] #16 #number of kpoints in the grid
OpticalIndex            = Materials[material]['OpticalIndex']
#print 

try: 
    kpoints_number_inp   = LookUpStringInFile("info", "Dimensions of the k-point grid")
    for i in np.arange(0,10):
        kpoints_number_inp=kpoints_number_inp.replace("  ", " ")
        if(verbose==1):
            print Header+"** Info: decreasing number of spaces"
            print kpoints_number_inp
    kpoints_number_array = kpoints_number_inp.split("=")
    kpoints_number_array = kpoints_number_array[1].split(" ")
    kpoints       = int(kpoints_number_array[1])
    print kpoints
except:
    kpoints                = kpoints_manual #manual input of k-point grid unit (assumes same in 3 directions)
    print Header+"** Error: failed to capture the number of k points automatically."
    #exit()

if(LogPlotting): 
    Efield_SI_log = np.linspace(np.log10(Efield_SI_min), np.log10(Efield_SI_max), num_fields)
    Efield_SI     = np.power(10., Efield_SI_log)
else:
    Efield_SI     = np.linspace(Efield_SI_min, Efield_SI_max, num_fields)


unitx=1E-4

for wavelength_SI in wavelength_SI_t:
    RefractiveIndex = OpticalIndex[wavelength_SI]
## Plot the splitted energy levels at Gamma.
    plt.figure()
    for MPI_number in MPI_numbers:
        for field in Efield_SI:
            eigenvalues, replicas = ComputeFloquetBandStructure(filename, nb_atoms, Z_electrons, kpoints, unocc_states, MPI_number, field, wavelength_SI, num_time_steps, listoffiles, Complex_valued_ME, field_polarization_dir, spin_occupation, verbose)
            # Plotting using external intensity, whence internal intensity acts on the sample
            field  = 0.5*c*epsilon_0*field**2*RefractiveIndex
            plt.scatter(np.ones(np.size(replicas))*field*unitx, Energy_Hartree_to_eV(np.real(replicas)), c="red", s=SizeOfDots, label="Floquet")
            
        Intensity_min = 0.5*c*epsilon_0*Efield_SI[0]**2*RefractiveIndex
        
        # Plotting eigenvalues at rest
        plt.scatter(np.ones(np.size(eigenvalues))*Intensity_min*unitx, Energy_Hartree_to_eV(np.real(eigenvalues)), c="black", s=10, label="GS")
        
        # Plotting dressed levels directly on the figure
        plt.scatter(np.ones(np.size(eigenvalues))*Intensity_min*unitx, Energy_Hartree_to_eV(np.real(eigenvalues))+MPI_number*h*c/wavelength_SI/e, c="red", s=10, label="replicates")
        plt.scatter(np.ones(np.size(eigenvalues))*Intensity_min*unitx, Energy_Hartree_to_eV(np.real(eigenvalues))-MPI_number*h*c/wavelength_SI/e, c="green", s=10, label="replicates")
        
        if(LogPlotting):
            plt.xscale('log')
            
        plt.title(r"$n=$"+str(round(MPI_number,0))+r", $\lambda=$"+str(int(wavelength_SI*1E9))+" nm, k="+str(kpoints)+r"$^3$")
        figurename="FloquetBandStructure-"+str(round(wavelength_SI*1E9,0))+"nm-k"+str(kpoints)+"MPI-"+str(round(MPI_number,0))
        #plt.xlabel("Field (V/m)")
        plt.xlabel(r"Laser intensity (W/cm$^{2}$)")
        plt.ylabel(r"Energy levels (eV)")
        plt.ylim(Materials[material]['ymin_ymax']) #Si
        #plt.ylim((-8,8)) #SiO2
        plt.tight_layout()
        plt.savefig(figurename+".eps")
        plt.savefig(figurename+".png")
        #plt.legend()
        plt.tight_layout()
    #plt.show()
