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
from libStark            import Kronecker, ComputeFloquetBandStructure

#Filenames = "ks_me_multipoles.*"
Filenames = "ks_me_dipole.*"
LocalFolder="./"
filename      = Filenames.split('.')[0]
print Filenames
# 1. Capture the static/dipole files. 
listoffiles=GetOrderedListOfFiles(Filenames)
print listoffiles

# 2. Check consistency of the number of files with k-point employed in the computation
print "Checking consistency on the number of files."
nb_atoms      = 1 #number of atoms #we actually treat a bulk, so it's one.
Z_electrons   = 4 #number of modeled valence electrons in the simulation box (3s)²(3p)² for Si atoms. 
kpoints                = 8 #k-point grid unit (assumes same in 3 directions)
unocc_states           = 4 #this is defined in the inp. If not, then it's zero. 
MPI_number             = 1
wavelength_SI          = 3200e-9
num_time_steps         = 100
num_fields             = 300 #number of electric field values for the plotting
Efield_SI_min          = 1E9
Efield_SI_max          = 1E11
SizeOfDots             = 0.5   #purely cosmetic!
Complex_valued_ME      = False
field_polarization_dir = 0 #0: x, 1: y, 2: z
spin_occupation        = 1 #1: levels with same spin should split, 2: levels of different spins dont split.

Efield_SI_log = np.linspace(np.log10(Efield_SI_min), np.log10(Efield_SI_max), num_fields)
Efield_SI     = np.power(10., Efield_SI_log)

OpticalIndex={ #data from Palik! 
    483e-9: 4.4028, 484e-9: 4.3964,
    800e-9: 3.6924167231905, 1600e-9: 3.4826, 
    2200e-9: 3.4548, 2600e-9: 3.4457, 3200e-9: 3.43438589865999}
RefractiveIndex = OpticalIndex[wavelength_SI]


## Plot the splitted energy levels at Gamma.
plt.figure()
plt.title("Laser dressing of electrons")
unitx=1E-4
for field in Efield_SI:
    eigenvalues, replicas = ComputeFloquetBandStructure(filename, nb_atoms, Z_electrons, kpoints, unocc_states, MPI_number, field, wavelength_SI, num_time_steps, listoffiles, Complex_valued_ME, field_polarization_dir, spin_occupation)
    # Plotting using external intensity, whence internal intensity acts on the sample
    field  = 0.5*c*epsilon_0*field**2*RefractiveIndex
    plt.scatter(np.ones(np.size(replicas))*field*unitx, Energy_Hartree_to_eV(np.real(replicas)), c="red", s=SizeOfDots, label="Floquet")
Intensity_min = 0.5*c*epsilon_0*Efield_SI[0]**2*RefractiveIndex
plt.scatter(np.ones(np.size(eigenvalues))*Intensity_min*unitx, Energy_Hartree_to_eV(np.real(eigenvalues)), c="black", s=10, label="GS")
plt.xscale('log')
plt.tight_layout()
figurename="FloquetBandStructure-"+str(round(wavelength_SI*1E9,0))+"nm-k"+str(kpoints)+"MPI-"+str(round(MPI_number,0))
plt.savefig(figurename+".eps")
#plt.xlabel("Field (V/m)")
plt.xlabel(r"Laser intensity (W/cm$^{2}$)")
plt.ylabel(r"Energy levels (eV)")
#plt.legend()
plt.tight_layout()
plt.show()
