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
from libStark            import Kronecker

Filenames = "ks_me_multipoles.*"
LocalFolder="./"
filename      = Filenames.split('.')[0]

# 1. Capture the static/multipoles files. 
listoffiles=GetOrderedListOfFiles(Filenames)
print listoffiles

# 2. Check consistency of the number of files with k-point employed in the computation
print "Checking consistency on the number of files."
nb_atoms   = 2 #number of atoms in simulation box
Z_electrons= 2 #number of modeled valence electrons in the simulation box
kpoints    = 8
m_points   = 3 #m=-1, 0, +1: 1 file for each. 
unocc_states=4

#Degenerated_Electrons=True #Do we describe 2 degenerated electron using 1 state? 
#if(Degenerated_Electrons): 
    ##True --> 0.5. False: 1. 
    #degenerated_electrons_f=0.5
#else:
    #degenerated_electrons_f=1
    
expected_num_files = nb_atoms * Z_electrons * m_points * kpoints**3
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

matrix_side=Z_electrons*nb_atoms+unocc_states
matrixelements=np.multiply(np.zeros((matrix_side, matrix_side)), 0e0+0e0j)
print "Info: size of coupling matrix elements: "+str(matrix_side)+"x"+str(matrix_side)

for file_id in gamma_point_files: #
    print "Info: opening file "+file_id
    # remove the commas from file_id
    ReplaceStringInFile(",", " ", file_id, "fixed_"+file_id)
    # capture the splitted real imag terms
    matrixelements_splitted = np.loadtxt("fixed_"+file_id, skiprows=4)
    for i in np.arange(0,matrix_side):
        for j in np.arange(0,matrix_side):
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
H_GS=np.zeros((matrix_side, matrix_side))*1j
#print np.size(H_GS[0])
for i in np.arange(0,np.size(H_GS[0])):
    H_GS[i,i]=eigenvalues[i]
print "Info: matrix of eigenvalues for the selected k-point. Matrix shape is "+str(np.shape(H_GS))

# 6.2: Add the dipolar matrix elements to the eigenvalues

## Function that returns the light-perturbed GS Hamiltonian for a specific k-point
def H_perturb(t, omega_AU, Efield_AU, H_GS,  matrixelements):
    return np.add(H_GS, Efield_AU*np.cos(omega_AU*t)*matrixelements)

# 6.3: introducing Fourier transform, multiphotonic levels and temporal averaging over one period.

Efield_SI = 1E9 #a.u.
wavelength_SI = 800e-9
omega_SI  = 2.*np.pi*c / wavelength_SI #this is in SI...
MPI_number=2
num_time_steps=1000

Efield_AU = Field_SI_to_AU(Efield_SI)
omega_AU  = Energy_eV_to_Hartree(omega_SI*hbar/e)

tmin=0.; tmax=2.*pi/omega_AU
dt=(tmax-tmin)/num_time_steps

# Initialization for t=0
print "Info: H(t)=H_GS + dipolar_matrix * Efield(t)"
H0=H_perturb(tmin, omega_AU, Efield_AU, H_GS, matrixelements)
print "Info: shape of H(t): "+str(np.shape(H0))

## Temporal integration of H_Floquet^{m,n}
def H_Floquet_mn_func(tmin, tmax, dt, omega_AU, MPI_number, n, Efield_AU, H_GS, matrixelements): #{{{
    H0_init = H_perturb(tmin, omega_AU, Efield_AU, H_GS, matrixelements) #just to initialize
    H0_sum = np.zeros(np.shape(H0_init))*0j #initialization
    for t in np.arange(tmin,tmax,dt): #integral
        H0=H_perturb(t, omega_AU, Efield_AU, H_GS, matrixelements) #H(t)
        H0_t1=omega_AU/2./np.pi * np.multiply(np.exp(1j*(MPI_number-n)*omega_AU*t), H0)*dt #e( i(m-n) omega t ) H(t)
        H0_sum = np.add(H0_sum, H0_t1) #have to integrate this
        
    H0_sum = H0_sum + Kronecker(MPI_number,n)*MPI_number*omega_AU
    return H0_sum
#}}}

H_Floquet_mn=H_Floquet_mn_func(tmin, tmax, dt, omega_AU, MPI_number, MPI_number, Efield_AU, H_GS, matrixelements)
#print H_Floquet_mn 
print "Info: element (m,n) of Floquet matrix constructed with success."
print "Info: Now, one has to write the full matrix."


print "Info: shape of element [m,n] of the matrix we want to write: "+str(np.shape(H_Floquet_mn))
H_Floquet_shape = (MPI_number*2+1) * (Z_electrons+unocc_states)
print "Info: shape of the complete Floquet matrix"
print H_Floquet_shape ** 2
H_Floquet=np.array((H_Floquet_shape, H_Floquet_shape)) #2d-list of arrays is needed
#print np.shape(H_Floquet)
for m in np.arange(-MPI_number,MPI_number+1,1):
    for n in np.arange(-MPI_number,MPI_number+1,1):
        H_Floquet_mn=H_Floquet_mn_func(tmin, tmax, dt, omega_AU, m, n, Efield_AU, H_GS, matrixelements) #tensor of 4th order...
        print "m="+str(m)
        print "m+MPI_number="+str(m+MPI_number)
        print "n="+str(n)
        print "n+MPI_number="+str(n+MPI_number)
        H_Floquet[m+MPI_number,n+MPI_number]=H_Floquet_mn
#exit()



# 6.4: Diagonalize the Floquet matrix
#w, v = LA.eig(H_perturb)
#print v

## 7. Plot the splitted energy levels at Gamma.
