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

Filenames = "ks_me_multipoles.*"
LocalFolder="./"
filename      = Filenames.split('.')[0]

# 1. Capture the static/multipoles files. 
listoffiles=GetOrderedListOfFiles(Filenames)
print listoffiles

# 2. Check consistency of the number of files with k-point employed in the computation
print "Checking consistency on the number of files."
nb_atoms      = 2 #number of atoms in simulation box
Z_electrons   = 2 #number of modeled valence electrons in the simulation box
kpoints       = 8
unocc_states  = 4
MPI_number    = 1
wavelength_SI = 3200e-9
num_time_steps= 100
num_fields    = 200
Efield_SI_min = 1E9
Efield_SI_max = 1E12
SizeOfDots    = 1

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
    eigenvalues, replicas = ComputeFloquetBandStructure(filename, nb_atoms, Z_electrons, kpoints, unocc_states, MPI_number, field, wavelength_SI, num_time_steps, listoffiles)
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
