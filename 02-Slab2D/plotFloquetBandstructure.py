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

Filenames = "ks_me_multipoles.*"
LocalFolder="./"
filename      = Filenames.split('.')[0]

# 1. Capture the static/multipoles files. 
listoffiles=GetOrderedListOfFiles(Filenames)
print listoffiles

# 2. Check consistency of the number of files with k-point employed in the computation
print "Checking consistency on the number of files."
nb_atoms   =2 #number of atoms in simulation box
Z_electrons=4 #number of modeled valence electrons in the simulation box
kpoints=8
m_points = 3 #m=-1, 0, +1: 1 file for each. 
unocc_states=4

Degenerated_Electrons=True #Do we describe 2 degenerated electron using 1 state? 
if(Degenerated_Electrons): 
    #True --> 0.5. False: 1. 
    degenerated_electrons_f=0.5
else:
    degenerated_electrons_f=1
    
expected_num_files = nb_atoms * Z_electrons * degenerated_electrons_f * m_points * kpoints**3
if(len(listoffiles) != expected_num_files):
    print "Error: set of data may be incomplete. Expected number of files: "+str(expected_num_files)
else: 
    print "Info: number of files is correct."

# 3. Construct the big matrice of dipolar transition inside the Floquet Hamiltonian. 
# WARNING: how to put m=-1, m=0, m=+1 ? does this order matter? 

# 3.1: NOTE: we start with gamma point only. 
# Gamma point index may be indicated in "info" file, at the gap indication. 
print "Capturing the Gamma point..."

print "** WARNING: this part is experimental!"
## METHOD 1: trusting the choice of octopus
info_file = "info"
print "Is there a direct gap value?"
try: 
    kpoints_index_inp   = LookUpStringInFile(info_file, "Direct gap")
    kpoints_index_array = kpoints_index_inp.split(" ")
    #print kpoints_index_array
    kpoints_index       = int(kpoints_index_array[4]) #BUG: dangerous...
    gamma_point_files_index=[kpoints_index*m_points, kpoints_index*m_points+1, kpoints_index*m_points+2]
    print "** Info: found a direct gap."
except:
    print "** Warning: could not find a direct gap. Looking for an indirect gap then."
    kpoints_index_inp   = LookUpStringInFile(info_file, "Indirect gap")
    kpoints_index_array = kpoints_index_inp.split(" ")
    kpoints_index       = int(kpoints_index_array[7]) #BUG: dangerous...
    gamma_point_files_index=[kpoints_index*m_points, kpoints_index*m_points+1, kpoints_index*m_points+2]

gamma_point_files = []
for num in gamma_point_files_index:
    gamma_point_file = filename+'.'+str(num)
    gamma_point_files.append(gamma_point_file)    
print "Files containing the matrix at Gamma point: "+str(gamma_point_files)

# 3.2. Capture the matrix in each file
# Structure is with complex numbers separated by spaces. 
# Commas are separating real and imag parts!
matrixelements=np.multiply(np.zeros((Z_electrons*nb_atoms, Z_electrons*nb_atoms)), 0e0+0e0j)

for file_id in gamma_point_files: #
    print "Info: opening file "+file_id
    # remove the commas from file_id
    ReplaceStringInFile(",", " ", file_id, "fixed_"+file_id)
    # capture the splitted real imag terms
    matrixelements_splitted = np.loadtxt("fixed_"+file_id, skiprows=4)
    for i in np.arange(0,Z_electrons*nb_atoms):
        for j in np.arange(0,Z_electrons*nb_atoms):
            matrixelements[i,j]=matrixelements_splitted[i, j]+1.j*matrixelements_splitted[i,j+1]
# 4. Plotting the dipolar transition matrices
    #plt.matshow(np.abs(matrixelements))
    #plt.title(file_id+", norm")
    #plt.colorbar()
    #plt.matshow(np.real(matrixelements))
    #plt.title(file_id+", real")
    #plt.colorbar()
    #plt.matshow(np.imag(matrixelements))
    #plt.title(file_id+", imag")
    #plt.colorbar()
#plt.show()
#print matrixelements
print "Info: successfull import of complex-valued dipolar matrix elements for m=-1, m=0, m=1."
print "** Warning: we still don't know how to combine these."
print "** Warning: therefore one selects only one file among m=-1,0,1 for now!"

## 5. Extract eigen values at k-point <kpoints_index> from "info" file
eigenvalues=GetEigenValue(kpoints_index, "info")
#print eigenvalues

## 6. Build then diagonalize the Floquet Hamiltonian. 
# 6.1:Build the eigenvalued matrix
H_GS=np.zeros((Z_electrons+unocc_states, Z_electrons+unocc_states))*1j
#print np.size(H_GS[0])
for i in np.arange(0,np.size(H_GS[0])):
    H_GS[i,i]=eigenvalues[i]
print H_GS

# 6.2: Add the dipolar matrix elements to the eigenvalues

Afield = 1 #a.u.
wavelength = 800e-9
t      = 0. #how to manage that? we have to integrate on time
omega  = 2.*np.pi*c / wavelength #this is in SI...

H_perturb=np.add(H_GS, Afield*np.cos(omega*t)*matrixelements)

plt.matshow(np.abs(H_perturb))
plt.colorbar()
plt.tight_layout()
plt.title(r"$H_{\mathrm{perturb}}$, norm")
plt.show()

# 6.3: introducing Fourier transform, multiphotonic levels and temporal averaging over one period.


# 6.4: Diagonalize the Floquet matrix
w, v = LA.eig(H_perturb)
print v

## 7. Plot the splitted energy levels at Gamma.
