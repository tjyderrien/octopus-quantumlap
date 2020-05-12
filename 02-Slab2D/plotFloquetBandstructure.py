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
from scipy.special import jn_zeros as BesselJZeros
from matplotlib.cm import get_cmap

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
material="Si" #SiO2 | Si

#exit()
wavelength_SI_t        = [3200e-9, 1600e-9, 800e-9, 484e-9, 400e-9] #484e-9, 343e-9] #10e-6

num_time_steps         = 100
num_fields             = 300 #number of electric field values for the plotting

Efield_SI_min          = IntensityToField(1E10*1E4)
Efield_SI_max          = IntensityToField(1E15*1E4)
Efield_AU_min          = Field_SI_to_AU(Efield_SI_min)
Efield_AU_max          = Field_SI_to_AU(Efield_SI_max)

Complex_valued_ME      = False
field_polarization_dir = 0 #0: x, 1: y, 2: z
spin_occupation        = 2 #1: levels with same spin should split, 2: levels of different spins don't split.

Renormalize            = False #True: limits the plots maximum as function of wavelength
PlotIntensity          = False  #True: plot with I. False: plot with E. 
SizeOfDots             = 0.5   #purely cosmetic!
LogPlotting            = 1
verbose                = 0 #Show info and debug messages. 
SingleKpoint           = 1 #1: plotting the whole band structure for a certain field

BesselRootsID_max      = 10


########### DONT TOUCH BELOW, thanks! #############

OpticalIndex_Si={ #NOTE: these are experimental data for Si (from Palik). They should be from linear TDDFT response !
    343e-9: 5.2749, 400e-9: 5.5674, 483e-9: 4.4028, 484e-9: 4.3964,
    800e-9: 3.6924167231905, 1600e-9: 3.4826, 
    2200e-9: 3.4548, 2600e-9: 3.4457, 3200e-9: 3.43438589865999, 10E-6: 3.4181}

OpticalIndex_SiO2={ #NOTE: these are experimental data for SiO2 (from refractiveindex.info). They should be from linear TDDFT response !
    343e-9: 1.5669, 400e-9: 1.5577, 483e-9: 1.5499, 484e-9: 1.5498,
    800e-9: 1.5383, 1600e-9: 1.5270, 
    2200e-9: 1.5175, 2600e-9: 1.5094, 3200e-9: 1.4996}

OpticalIndex_Mo={3200e-9: 1.6626, 1600e-9: 1.5686, 1030e-9: 2.4201, 800e-9: 3.6371, 400e-9: 0.42388}

Si={'MPI_numbers': [1,2,3], 'kpoints_manual': 16, 
    'nb_atoms': 1, 'unocc_states': 4, 
    'Z_electrons': 4, 'OpticalIndex': OpticalIndex_Si, 
    'ymin_ymax': (3,8), 'CaptureDirectGap': True}

SiO2={'MPI_numbers': [1,2], 'kpoints_manual': 16, 
    'nb_atoms': 1, 'unocc_states': 6, 
    'Z_electrons': 24, 'OpticalIndex': OpticalIndex_SiO2, 
    'ymin_ymax': (-8,16), 'CaptureDirectGap': True}

Mo={'MPI_numbers': [1], 'kpoints_manual': 16, 
    'nb_atoms': 1, 'unocc_states': 4, 
    'Z_electrons': 11-4, 'OpticalIndex': OpticalIndex_Mo, 
    'ymin_ymax': (-30,30), 'CaptureDirectGap': False}

Materials={'Si': Si, 'SiO2': SiO2, 'Mo': Mo} #this contains all materials properties we need to automatize

# 2. Check consistency of the number of files with k-point employed in the computation
nb_atoms                = Materials[material]['nb_atoms'] #1 #number of atoms #we actually treat a bulk, so it's one.
Z_electrons             = Materials[material]['Z_electrons'] #SiO2: 24 electrons. #Si: 4 #number of modeled valence electrons in the simulation box 
# (3p)² x 2 for Si (biatomic) lattice. 
# SiO2 ? matrix is 30 e- states, 6 unocc. 
unocc_states            = Materials[material]['unocc_states'] #4 #SiO2: 6. #Si: 4. This is defined in the inp. If not, then it's zero. 

CaptureDirectGap        = Materials[material]['CaptureDirectGap']

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


print "Info: Capturing the file for the transition of interest."
gamma_point_files, gamma_point_files_index = ExtractFileForGivenKpoint(filename, kpoints, listoffiles, CaptureDirectGap, verbose, field_polarization_dir)

matrixelements, matrix_side = ReadDipolarMatrixElements(Z_electrons, unocc_states, gamma_point_files, Complex_valued_ME)

## 5. Extract eigen values at k-point <kpoints_index> from "info" file
eigenvalues=GetEigenValue(gamma_point_files_index[0], "info", verbose)
if(verbose==1):
    print eigenvalues

#tunneling_stops=[]

for wavelength_SI in wavelength_SI_t:
    
    RefractiveIndex = OpticalIndex[wavelength_SI]
    
    ## Renormalize the intensity range using refractive index
    if(Renormalize): 
        Efield_SI_max = Efield_SI_max/np.sqrt(OpticalIndex[wavelength_SI])
    
    if(LogPlotting): 
        Efield_SI_log = np.linspace(np.log10(Efield_SI_min), np.log10(Efield_SI_max), num_fields)
        Efield_SI     = np.power(10., Efield_SI_log)
    else:
        Efield_SI     = np.linspace(Efield_SI_min, Efield_SI_max, num_fields)

## Plot the splitted energy levels at Gamma.
    for MPI_number in MPI_numbers:
        plt.figure()
        for field in Efield_SI:
            eigenvalues, replicas, dipoles = ComputeFloquetBandStructure(filename, nb_atoms, Z_electrons, kpoints, unocc_states, MPI_number, field, wavelength_SI, num_time_steps, eigenvalues, matrixelements, Complex_valued_ME, field_polarization_dir, spin_occupation, verbose)
            # Plotting using external intensity, whence internal intensity acts on the sample
            if(PlotIntensity):
                unitx=1E-4
                field  = 0.5*c*epsilon_0*field**2*RefractiveIndex
            else:
                unitx=1E0
                field = field*np.sqrt(RefractiveIndex)
            plt.scatter(np.ones(np.size(replicas))*field*unitx, Energy_Hartree_to_eV(np.real(replicas)), c="red", s=SizeOfDots, label="Floquet replicates") #these are the interesting lines for the final plot
            
            # TODO: Prepare the table [wavelength, n, E, Rabi Freq such as BesselJ(n,RabiFreqs_SI/omegaLaser)=0. #BruteForce method. 
            # TODO: Do we need the function BesselJZeros(n,x)? To search omegaLaser that induce a stop in tunneling! 
            #tunneling_stops.append(np.sqrt(dipoles * field_AU / BesselJZeros(n,k)))
            
        if(PlotIntensity):
            Intensity_min = 0.5*c*epsilon_0*Efield_SI[0]**2*RefractiveIndex
        else: 
            Intensity_min = Efield_SI[0]*np.sqrt(RefractiveIndex)
        
        # Plotting eigenvalues at rest
        plt.scatter(np.ones(np.size(eigenvalues))*Intensity_min*unitx, Energy_Hartree_to_eV(np.real(eigenvalues)), c="black", s=10, label="GS eigenvalues")
        
        if(LogPlotting):
            plt.xscale('log')
            
        plt.title(r"$n=$"+str(round(MPI_number,0))+r", $\lambda=$"+str(int(wavelength_SI*1E9))+" nm, k="+str(kpoints)+r"$^3$")
        figurename="FloquetBandStructure-"+str(round(wavelength_SI*1E9,0))+"nm-k"+str(kpoints)+"MPI-"+str(round(MPI_number,0))
        if(not LogPlotting):
            figurename+="NoLog"
        if(not PlotIntensity): 
            figurename+="-Field"
            plt.xlabel("Field (V/m)")
        else:
            plt.xlabel(r"Laser intensity (W/cm$^{2}$)")
        plt.ylabel(r"Energy levels (eV)")
        plt.ylim(Materials[material]['ymin_ymax']) #Si
        #plt.ylim((-8,8)) #SiO2
        plt.tight_layout()
        plt.savefig(figurename+".eps")
        plt.savefig(figurename+".png")
        # Plotting dressed levels directly on the figure, from n=1 to MPI_number
        for MPI in np.arange(1,MPI_number+1):
            plt.scatter(np.ones(np.size(eigenvalues))*Intensity_min*unitx, Energy_Hartree_to_eV(np.real(eigenvalues))+MPI*h*c/wavelength_SI/e, c="red", s=10, label="non-interacting replicates")
            plt.scatter(np.ones(np.size(eigenvalues))*Intensity_min*unitx, Energy_Hartree_to_eV(np.real(eigenvalues))-MPI*h*c/wavelength_SI/e, c="blue", s=10, label="non-interacting replicates")
        
        plt.savefig(figurename+"-Suppl"+".eps")
        plt.savefig(figurename+"-Suppl"+".png")
        #plt.legend()
    #plt.show()
#print Header+"Preparing (omega,Efield) values that should disable the tunneling in the sample..."
#plt.figure()

print "=== Assessing tunneling stopping possibilities... ==="
Efield_AU_log = np.linspace(np.log10(Efield_AU_min), np.log10(Efield_AU_max), num_fields)
Efield_AU     = np.power(10., Efield_AU_log)

for MPI_number in np.arange(-np.max(MPI_numbers), np.max(MPI_numbers)+1): #[Kato, PRB 100 081203(R) tells to use -MPI_numbers..MPI_numbers]
    tunneling_stops=[]
    
    figurename = str(material)+"-MatrixElements-k-"+str(kpoints)+"MPI-"+str(round(MPI_number,0))
    #for BesselRootsID in np.arange(1,BesselRootsID_max):
    BesselZeros=BesselJZeros(MPI_number,BesselRootsID_max)
    #print "BesselRoot #"+str(BesselRootsID_max)
    print "BesselJZeros("+str(MPI_number)+", "+str(BesselRootsID_max)+")="+ str(BesselZeros)
    for RootNumber in np.arange(0,np.size(BesselZeros)): #NOTE: we mix the various roots of Bessel function
        tunneling_stops.append(np.divide(matrixelements, BesselZeros[RootNumber])) # = omega**2/E. #NOTE: concatenation of all the roots of Bessel. 
    print np.shape(tunneling_stops)
    #print tunneling_stops #omega**2/E (in at.u.)
    #plt.figure()
    plt.matshow(np.log10(tunneling_stops[0])) #omega**2/E
    plt.title(r'$n=$'+str(MPI_number)) #+", slice "+str(BesselRootsID_max-1))
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(figurename+".eps")
    plt.savefig(figurename+".png")
    
    figurename=str(material)+"-BesselTunnelingStop-"+str(kpoints)+"MPI-"+str(round(MPI_number,0))+"-BesselRoot-"+str(BesselRootsID_max)
    
    print "Plotting couples of (E, omega) that blocks transition with "+str(MPI_number)+" photons."
    #for wavelength_SI in wavelength_SI_t: #for each wavelength
    plt.figure()
    
    # Number of colors we need. 
    NumberOfElements=(matrix_side/2)
    cm = get_cmap('hsv')
    colors=[cm(1.*i/(NumberOfElements)) for i in range(NumberOfElements)]
    xy = range(NumberOfElements)
    colorlist=[colors[x] for x in xy]
    
    plot1=[]
    for trans_i in np.arange(0,matrix_side/2): #transition from i
        for trans_j in np.arange(matrix_side/2,matrix_side,1): #transition to j
            for RootNumber in np.arange(0,np.size(BesselZeros)):
                plot_field = plt.scatter(Field_AU_to_SI(Efield_AU), Energy_Hartree_to_eV(np.sqrt(tunneling_stops[RootNumber][trans_i][trans_j]*Efield_AU)), s=10*SizeOfDots, label=str(trans_i)+r"$\rightarrow$"+str(trans_j), c=colorlist[trans_j-matrix_side/2]) #colorlist[trans_j] is also possible #*matrix_side+trans_j]) #energies [eV]
            plot1+=[plot_field] #stores the last plot caption in a series of roots
    
    labels1 = [l.get_label() for l in plot1]
    plt.ylabel(r"$\hbar \omega$ (eV)")
    #plt.ylabel(r"$\lambda$ (m)")
    plt.legend(plot1, labels1, loc="best", ncol=2)
    plt.xlabel(r"$E$ (V/m)")
    plt.xscale('log')
    plt.ylim((0,4))
    plt.title(r"Tunneling $(E,\omega)$, $n=$"+str(MPI_number))
    plt.savefig(figurename+".eps")
    plt.savefig(figurename+".png")
    #plt.show()
