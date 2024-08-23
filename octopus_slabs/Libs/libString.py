#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2013-2024 T. J.-Y. Derrien
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

## @package libString
## Look for strings in files.

import glob
# from ase.io.octopus import read_static_info_kpoints
# IMPORT LIBRARIES
import os
# from libDatabase import * #uses symbolic link
# from ../../../spp_extended_theory/SimpleSPP/libDatabase import FilterDatabaseContains
# from importlib import import_module

# libDatabase = import_module("spp_extended_theory.SimpleSPP.libDatabase")
from spp_extended_theory.Libs.libDatabase import *
# from octopus_slabs.Libs.libPlotting import CaptureKPointGridSizeFromINP


## Returns the name of the most recent octopus logfile
def LookUpLatestLogFile(Filename): #{{{{
  Header = "[libString: LookUpLatestLogFile] "
  import glob
  newest = max(glob.iglob(Filename+'*.log'), key=os.path.getctime)
  print(Header+"List of files: "+str(newest))
  # We shall split the path by slashes. 
  newest_splitted = newest.split("/")
  print(Header+"List of files: "+str(newest_splitted))
  newest_file = newest_splitted[len(newest_splitted)-1] #last element of the path is selected
  print(Header+"List of files: "+str(newest_splitted))
  result = [int(s) for s in newest_file.split(".") if s.isdigit()]
  print(Header+"List of files: "+str(result))
  # For some folders, result was empty (?!)
  try:
    final = max(result)
  except: 
    final = result
  return final
#}}}

## Check whether a file contains binary data
def contains_binary_data(file_path):
  with open(file_path, 'rb') as file:
      while True:
          byte = file.read(1)
          if not byte:
              break
          if byte < b' ' and byte not in (b'\t', b'\n', b'\r'):
              return True
  return False

## Looks for a string in a file. 
# Returns the full line where the element is found. 
# If duplicate files are found, the most recent file is selected.
# Stops as soon as first element is found, to avoid taking next lines
def LookUpStringInFile(Filename, string): #{{{
  import glob
  newest = max(glob.iglob(Filename), key=os.path.getctime)

  #Open the log file
  searchfile = open(newest, "r")

  # Search into this file the <line> with <string>
  result=''

  for line in searchfile:
    if string in line:
      result = line;
      break; #stop as soon as first element is found, to avoid taking next lines
  return result
#}}}

## Return the line that follows an identified pattern
def LookUpStringReturnPreviousLine(Filename, string):  # {{{
    import glob
    newest = max(glob.iglob(Filename), key=os.path.getctime)

    # Open the log file
    searchfile = open(newest, "r")

    # Search into this file the <line> with <string>
    result = ''
    found = False
    linePrev = ""
    for line in searchfile:
        if string in line:
            found = True
        if (found):
            result = linePrev  # returns the line N+1 when pattern is found in line N
            break
        linePrev = line
    return result
    # }}}

## Return the line that follows an identified pattern
def LookUpStringReturnNextLine(Filename, string): #{{{
  import glob
  newest = max(glob.iglob(Filename), key=os.path.getctime)
  
  #Open the log file
  searchfile = open(newest, "r")

  # Search into this file the <line> with <string>
  result=''
  found=False
  for line in searchfile:
      if(found):
          result=line #returns the line N+1 when pattern is found in line N
          break
      if string in line:
          found=True
  return result
  #}}}

## Looks for subfolders 
def FindAllDirectories(directory): 
  for dirname, dirnames, filenames in os.walk(directory):
         # print path to all subdirectories first.
         for subdirname in dirnames:
                print((os.path.join(dirname, subdirname)))

         ## print path to all filenames.
         #for filename in filenames:
                #print(os.path.join(dirname, filename))

         # Advanced usage:
         # editing the 'dirnames' list will stop os.walk() from recursing into there.
         if '.git' in dirnames:
                # don't go into any .git directories.
                dirnames.remove('.git')
  WholeTable=dirname+subdirname
  return WholeTable

## Generate subfolder list and filters using the provided KeyWord. 
# Equivalent of <find . | grep 'keyword'> in Bash. 
def FilterDirectories(keyword, path, FieldIndex):
  ListOfFolders = FindAllDirectories(path)
  ListOfFoldersFiltered = FilterDatabaseContains(ListOfFolders, keyword, FieldIndex)
  return ListOfFoldersFiltered

## Returns the last line of file
# Returns the last line (<string>) of a file
# Source: http://stackoverflow.com/questions/3346430/what-is-the-most-efficient-way-to-get-first-and-last-line-of-a-text-file
# @param file <string> filename
def GetLastLine(file):
  with open(file, "rb") as f:
         first = f.readline()      # Read the first line.
         f.seek(-2, 2)             # Jump to the second last byte.
         while f.read(1) != b"\n": # Until EOL is found...
                f.seek(-2, 1)         # ...jump back the read byte plus one more.
         last = f.readline()       # Read last line.
         return last

## Splitting, sorting and re-packing the table of ks_me_dipoles from octopus.
def GetOrderedListOfFiles(Filename): #{{{
    filelist = glob.iglob(Filename) #makes selection over files, but not sorted
    #NOTE: names are like ks_me_dipole.k999_y
    splitted_t = [[],[], []]
    names=[]; numbers=[];  directions=[]
    for file in filelist: #{{{
        splited=file.split(".")
        names.append(splited[0])
        splited_2 = (splited[1].split("_"))
        directions.append(splited_2[1])
        splited_3 = splited_2[0].split("k")
        numbers.append(int(splited_3[1]))
    #}}}
    # Sorting the table properly
    lists = sorted(zip(*[names, numbers, directions]))
    names_sorted, numbers_sorted, directions_sorted = list(zip(*lists))
    #numbers.sort()
    filelist = np.transpose([names_sorted, numbers_sorted, directions_sorted])

    finallist = []
    for line in filelist:
        newline=str(line[0])+".k"+str(line[1])+"_"+str(line[2])
        finallist.append(newline)
    return finallist
#}}}


## Replaces a string from a source file in a new file
# @param string: pattern to be replaced
# @param string: new string replacing the pattern
# @param file_inp: name of source file
# @param file_out: name of target file
def ReplaceStringInFile(string, replace_with, file_inp, file_out):
    Header="libString: ReplaceStringInFile(): "
    #print Header+"input file: "+file_inp
    with open(file_inp, "rt") as fin:
        with open(file_out, "wt") as fout:
            for line in fin:
                fout.write(line.replace(string, replace_with))
    return 0

## Gets Fermi energy and returns its value in Hartree
def GetFermiLevel(filename):
    Capture=LookUpStringInFile(filename, "Fermi energy")
    FermiEnergy=float(Capture.split("=")[1].strip().split(" ")[0])
    return FermiEnergy

## Parser by block
# @param some_source: source file name
# @param header: pattern to find the header
def group_by_heading( some_source, header ):
    buffer= []
    for line in some_source:
        if line.startswith( header ):
            if buffer: yield buffer
            buffer= [ line ]
        else:
            buffer.append( line )
    yield buffer
    

##Getting k-grid value: capture the first line after the declaration on k #NOTE: this may be already done in ASE
    # Warning: This assumes isotropic k-point grid meshing. Can be improved by necessity in future.
def CaptureKPointGridSizeFromINP(filename):
    kpointgrid_inp=LookUpStringReturnNextLine(filename, "%KPointsGrid")
    kpointgrid_inp2=kpointgrid_inp.strip()
    kpointgrid_inp3=kpointgrid_inp2.split("|")
    kpointgrid_inp4=kpointgrid_inp3[0] #capture the first coefficient of the grid
    kpointgrid=float(kpointgrid_inp4)
    del kpointgrid_inp, kpointgrid_inp2, kpointgrid_inp3, kpointgrid_inp4
    return kpointgrid

## Capture the k-point coordinate (between 0 and 1) for plotting
def CapturePlottingCoordinateOfKpointsPathPositions(filename="bandstructure"):
    epsilon = 7E-1
    import numpy as np
    import matplotlib.pyplot as plt
    bands = np.loadtxt(filename, dtype="float32", skiprows=1)
    kindex = bands[:, 0]
    kx = bands[:, 1]
    ky = bands[:, 2]
    kz = bands[:, 3]
    slope_x = np.gradient(kx) / np.gradient(kindex)
    slope_y = np.gradient(ky) / np.gradient(kindex)
    slope_z = np.gradient(kz) / np.gradient(kindex)
    discont_x = np.fabs(np.gradient(slope_x) / np.gradient(kindex))
    discont_y = np.fabs(np.gradient(slope_y) / np.gradient(kindex))
    discont_z = np.fabs(np.gradient(slope_z) / np.gradient(kindex))
    condition = np.logical_or(np.fabs(discont_x) > epsilon, np.fabs(discont_y) > epsilon)
    condition2 = np.logical_or(condition, np.fabs(discont_z) > epsilon)
    test2 = kindex[(condition2)]
    kindex_dis = np.reshape(test2, [int(len(test2) / 3), 3])
    result = np.mean(kindex_dis, axis=1)
    return(result)

## Get the k-grid size from the second line containing Gamma point
# @param filename: path/filename of the bandstructure file
# Assumes that KPointPath starts with Gamma point
# Assumes that k-grid size is same in every dimensions.
# WRONGLY assumes that bandstructure follows k_idx. It's not the case.
def CaptureKPointGridSizeFromBandStructureFile(filename="bandstructure"):
    bands = np.loadtxt(filename, skiprows=1)
    kpointpath_index=bands[:,0]
    kpointpath_kx   =bands[:,1]
    kpointpath_ky   =bands[:,2]
    kpointpath_kz   =bands[:,3]
    Mask=np.logical_and(np.logical_and(np.logical_and(kpointpath_kx==0, kpointpath_ky==0), kpointpath_kz==0), kpointpath_index==0)
    print(np.where(Mask))
    return np.where(Mask) #returns the line number of last line before gamma point
    #kgrid=LookUpStringReturnPreviousLine(filename, "0,0,0")

## Captures a block of instructions
# TODO: This function should be replaced by a parser for reading octopus file (ASE must have this...)
# NOTE: yes, and the function in ase.io.octopus is named read_static_info_kpoints()
def CaptureBlock(filename, string, shift=0):
    Header="libString: CaptureBlock(): "
    import pathlib
    linenumber=0; loopshift=0
    block=[]
    CaptureText=False
    filepath=pathlib.Path(filename)
    with filepath.open("r") as file_reader:
        for line in file_reader:
            # print(Header+"** Debug: "+line)
            linenumber+=1
            if(string in line): #if we find our pattern, we capture the text that follows
                CaptureText = True
#                print(Header+"** Info: text found line "+str(linenumber)+". Acquisition can start. Loopshift = "+str(loopshift))
            if("\n" == line and CaptureText):
                CaptureText = False
#                print(Header+"** Info: acquisition can stop.")
            if(CaptureText):
                if(loopshift<shift):
                    loopshift += 1
                else:
                    block.append(line.strip().split()) #building up the block
            #Then we can capture the block until finding two newlines finishing the block
    return(block)

def read_static_info_kpoints(fd):
    for line in fd:
        if line.startswith('List of k-points'):
            break

    tokens = next(fd).split()
    assert tokens == ['ik', 'k_x', 'k_y', 'k_z', 'Weight']
    bar = next(fd)
    assert bar.startswith('---')

    kpts = []
    weights = []

    for line in fd:
        # Format:        index   kx      ky      kz     weight
        m = re.match(r'\s*\d+\s*(\S+)\s*(\S+)\s*(\S+)\s*(\S+)', line)
        if m is None:
            break
        kxyz = m.group(1, 2, 3)
        weight = m.group(4)
        kpts.append(kxyz)
        weights.append(weight)

    ibz_k_points = np.array(kpts, float)
    k_point_weights = np.array(weights, float)
    return dict(ibz_k_points=ibz_k_points, k_point_weights=k_point_weights)

## Get the k-grid size from the second line containing Gamma point
# @param filename: path/filename of the bandstructure file
# Assumes that KPointPath starts with Gamma point
# Assumes that k-grid size is same in every dimensions.
def CaptureKPointGridSizeFromInfoFile(filename="info"):
    shift=3 #3 lines between signal and text to capture
    block=CaptureBlock(filename, "List of k-points:", shift) #line of first k-point
    block_t = np.asfarray(block)
    kpointpath_index=block_t[:,0]
    kpointpath_kx   =block_t[:,1]
    kpointpath_ky   =block_t[:,2]
    kpointpath_kz   =block_t[:,3]
    occupations     =block_t[:, 4]

    """"
    # WARNING!!! Hypothesis employed here is that k-point path starts with a Gamma point.
    test1 = (kpointpath_kx == 0)
    test2 = (kpointpath_ky == 0)
    test3 = (kpointpath_kz == 0)
    Mask = np.logical_and(np.logical_and(test1, test2), test3)
    """

    test4 = (occupations == 0) #this captures the first point of the k-point-path
    Mask = test4
    # not GammaPoint and occupations == 0 # we start the k-point-path (by something else than Gamma)
    # GammaPoint and occupations == 0 # we are in the k-point-path, but maybe not at the beginning (several Gammas are present in the path)
    # not GammaPoint and occupations == 1 # we are still inside 3D k-point-grid
    # GammaPoint and occupations == 1 # we are still inside the 3D k-point-grid
    #Mask = np.logical_and(np.logical_and(np.logical_and(kpointpath_kx==0, kpointpath_ky==0), kpointpath_kz==0), kpointpath_index==0)
    # block = read_static_info_kpoints(filename)

    NumberOf3dKpoints=np.where(Mask)[0][1]



    return round(np.power(NumberOf3dKpoints,1/3))

## Capture eigenvalues for a given k-point from "info" file
# @param kpoints_index: index of the k-point of interest
# @param NumberOfEigenValues: limits the capture of number of eigenvalues.
# @param filename: "info" file.
def GetEigenValue(kpoints_index, NumberOfEigenValues, filename="info", verbose=1):
     import re
     import numpy as np
     Header = "[libString: GetEigenValue(): ]"
     pattern = "#k ="
     if verbose == 1:
         print("Info: Construction of the eigen-values matrix.")
         print("Pattern of interest: \"" + pattern + "\"")

     eigenvalues = []

     with open(filename, "r") as source:
         for heading_and_lines in group_by_heading(source, pattern):
             heading = heading_and_lines[0]
             lines = heading_and_lines[1:]

             # Define the regular expression pattern
             prog = re.compile(r"#\s*k\s*=\s*(\d+)\s*,")

             # Search for the pattern in the heading
             result = prog.search(heading)

             if result:
                 k_value = int(result.group(1))  # Extract the kpoints_index

                 # Check if k_value matches the desired kpoints_index
                 if k_value == kpoints_index:
                     if verbose == 1:
                         print(kpoints_index)
                     break  # Stop iteration when pattern is found

     if verbose == 1:
         print(Header + str(heading))
         print(Header + str(lines[0:NumberOfEigenValues]))

     # Extract eigenvalues from the lines
     for line in lines[:NumberOfEigenValues]:
         data = line.split()
         if len(data) >= 3:
             eigenvalues.append(float(data[2]))  # Convert to float and append to the list

     # Convert the eigenvalues list to a NumPy array and return it
     return np.array(eigenvalues)

## Optimization from ChatGPT 4
import re
import numpy as np
# from joblib import Parallel, delayed
from tqdm import tqdm

def load_data(filename):
    pattern = "#k ="
    data_by_k = {}
    current_k = None
    current_lines = []

    with open(filename, "r") as source:
        for line in source:
            if pattern in line:
                if current_k is not None:
                    data_by_k[current_k] = current_lines
                match = re.search(r"#\s*k\s*=\s*(\d+)\s*,", line)
                if match:
                    current_k = int(match.group(1))
                    current_lines = [line]
            elif current_k is not None:
                current_lines.append(line)

        if current_k is not None:
            data_by_k[current_k] = current_lines

    return data_by_k


def GetEigenValueFromData(data, kpoints_index, NumberOfEigenValues, verbose=1):
    Header = "[libString: GetEigenValue(): ]"

    if verbose == 1:
        print("Info: Construction of the eigen-values matrix.")
        print(f"Looking for kpoints_index: {kpoints_index}")

    eigenvalues = []
    if kpoints_index in data:
        lines = data[kpoints_index]

        if verbose == 1:
            print(Header + str(lines[0]))
            print(Header + str(lines[1:NumberOfEigenValues + 1]))

        # Extract eigenvalues from the lines
        for line in lines[1:NumberOfEigenValues + 1]:  # Skip the header line
            line = line.strip()
            data_fields = line.split()
            if len(data_fields) >= 3:
                eigenvalues.append(float(data_fields[2]))  # Convert to float and append to the list

    return np.array(eigenvalues)


def process_kidx(data, k_index, kidx, NumberOfBands):
    return GetEigenValueFromData(data, int(k_index[int(kidx)]), NumberOfBands, verbose=0)


def GetEigenValues(filename, k_index, bandstructure_pd, NumberOfBands):
    data = load_data(filename)  # Load and structure data

    """ 
    # Use parallel processing with tqdm for progress bar
    
    with tqdm(total=len(bandstructure_pd.k_index.index)) as pbar:
        Eigenvalues_list = Parallel(n_jobs=-1)(
            delayed(lambda kidx: process_kidx(data, k_index, kidx, NumberOfBands) or pbar.update(1))
            (kidx) for kidx in bandstructure_pd.k_index.index
        )
    """

    Eigenvalues_list = []
    for kidx in tqdm(bandstructure_pd.k_index.index):
        eigenvalues = GetEigenValueFromData(data, int(k_index[int(kidx)]), NumberOfBands, verbose=0)
        Eigenvalues_list.append(eigenvalues)

    return Eigenvalues_list

## End of code optimized by ChatGPT

## Among a list of files, extract the single-direction dipolar-date filename that matches with: direct transition,
## gamma point xor for the given list of k-points (indicated by their index).
# @param filename: common part of the multiple filenames
# @param kpoints: indicative number of 3D kpoints in the simulation (not the path)
# @pram listoffiles: full list of the existing dipolar files
# @param listofkpoints: list of interesting kpoints
# @param CaptureDirectGap: bool for automatic detection.
def ExtractFileForGivenKpoints(filename, kpoints, listoffiles, listofkpoints, CaptureDirectGap=True, verbose=False, field_polarization_dir=0): #{{{
    Header="[libString: ExtractFileForGivenKpoints()] "                                                                                                                                                    
    spatial_directions      = 3 #1 file for each direction. 
    #field_polarization_dir  = 0 #0: x, 1: y, 2: z. 
    directions={0: "_x", 1: "_y", 2: "_z"}
    #Degenerated_Electrons=True #Do we describe 2 degenerated electron using 1 state? 
    #if(Degenerated_Electrons): 
        ##True --> 0.5. False: 1. 
        #degenerated_electrons_f=0.5
    #else:
        #degenerated_electrons_f=1
    
    # 2. Check consistency of the number of files with k
    #expected_num_files = spatial_directions * kpoints**3 #BUG: will never be good, since KpointPath is added after kPointGrid. Then we'll always fall into case of plotting all kpointpath. 
    #if(len(listoffiles) != expected_num_files):
        #print(Header+"Warning: set of data may be incomplete. Expected number of files: "+str(expected_num_files))
    #else: 
        #if(verbose==1): 
            #print(Header+"Info: number of files is correct.")
        #else: 
            #print("")

    # 3. Construct the big matrice of dipolar transition inside the Floquet Hamiltonian. 
    # WARNING: how to put m=-1, m=0, m=+1 ? does this order matter? 

    # 3.1: NOTE: we start with gamma point only. 
    # Gamma point index may be indicated in "info" file, at the gap indication. 
    
    #print "** WARNING: this part is experimental!"
    ## METHOD 1: trusting the choice of octopus
    info_file = "info"
    
    # NOTE: eigen values are sorted by norm. Then taking the first k point should be sufficient. 
    if(CaptureDirectGap):
        print(Header+"Capturing the Gamma point.")
        print(Header+"Is there a direct gap value? Useful only for band-gap materials.")
        try: #we try to capture direct gap from octopus static/info file.
            kpoints_index_inp   = LookUpStringInFile(info_file, "Direct gap")
            for i in np.arange(0,10):
                kpoints_index_inp=kpoints_index_inp.replace("  ", " ")
                if(verbose==1):
                    print(Header+"** Info: decreasing number of spaces")
                    print(Header+" "+str(kpoints_index_inp))
            kpoints_index_array = kpoints_index_inp.split(" ")
            #print kpoints_index_array
            kpoints_index       = int(kpoints_index_array[4])
            gamma_point_files_index=[kpoints_index]
            print(Header+"** Info: found a direct gap: "+str(gamma_point_files_index))
        except: #we select Gamma point by hand
            print(Header+"** WARNING: could not find a direct gap. Fetching the Gamma point.")
            print(Header+"Selecting the first k-point, corresponding to Gamma.")
            #kpoints_index_inp   = LookUpStringInFile(info_file, "k #1")
            #for i in np.arange(0,10):
                #kpoints_index_inp=kpoints_index_inp.replace("  ", " ")
                #if(verbose==1):
                    #print Header+"** Info: decreasing number of spaces"
                    #print kpoints_index_inp
            #kpoints_index_array = kpoints_index_inp.split(" ")
            ##print kpoints_index_array
            print(Header+"** Warning: kpoints_index selection is MANUAL.")
            kpoints_index       = 1
            gamma_point_files_index=[kpoints_index]
            #exit()
            #kpoints_index_inp   = LookUpStringInFile(info_file, "Indirect gap")
            #kpoints_index_array = kpoints_index_inp.split(" ")
            #kpoints_index       = int(kpoints_index_array[7]) #BUG: dangerous...
            #gamma_point_files_index=[kpoints_index]
            ##gamma_point_files_index=[kpoints_index*spatial_directions, kpoints_index*spatial_directions+1, kpoints_index*spatial_directions+2]
    else: #we will here output all k-points corresponding to the direction of polarization. 
        if(verbose==1):
            print((Header+"** Info: capturing all the indicated k-points..."))
            print((Header+"** Info: number of selected k points: "+str(len(listofkpoints))))
        gamma_point_files_index=listofkpoints
        
        #print("** Info: number of k points: "+str(len(listoffiles)/3.))
        #gamma_point_files_index=np.arange(1,len(listoffiles)/3+1) #NOTE: name is misleading. We don't want only gamma points here but all k-points.
        
    # Managing filenames for dipolar matrix elements
    gamma_point_files = []
    for num in gamma_point_files_index:
        gamma_point_file = filename+'.k'+str(num)+directions[field_polarization_dir]
        gamma_point_files.append(gamma_point_file)
    if(verbose==1):
        print(Header+"Files containing the matrix at desired point: "+str(gamma_point_files)) #not only for Gamma point
    return gamma_point_files, gamma_point_files_index
#}}}

# Capture the matrix in each file
# Structure is with real-valued numbers separated by spaces. 
# Commas are separating real and imag parts.
def ReadDipolarMatrixElements(Z_electrons, unocc_states, gamma_point_files, Complex_valued_ME=False, verbose=False, spin_occupation=2): #{{{
    Header="[libStark] ReadDipolarMatrixElements(): "
    matrix_side=Z_electrons+unocc_states
    matrixelements=np.multiply(np.zeros((matrix_side, matrix_side)), 0e0+0e0j)
    if(verbose==1): 
        print("Info: size of coupling matrix elements: "+str(matrix_side)+"x"+str(matrix_side))

    for file_id in gamma_point_files: #gamma_point_files should be an array of strings
        print(("** Info: opening file "+file_id))
        if(verbose==1): 
            print("Info: opening file "+file_id)
        if(Complex_valued_ME): 
            # remove the commas from file_id
            # NOTE: This is in case of complex-valued matrix from ks_multipoles
            ReplaceStringInFile(",", " ", file_id, "fixed_"+file_id)
            # capture the splitted real imag terms
            matrixelements_splitted = np.loadtxt("fixed_"+file_id, skiprows=4)
            for i in np.arange(0,matrix_side):
                for j in np.arange(0,matrix_side):
                    matrixelements[i,j]=matrixelements_splitted[i, j]+1.j*matrixelements_splitted[i,j+1]
        else: #case for the ks_me_dipole
            if(verbose==1):
                print(Header+"Info: opening "+file_id)
            matrixelements=np.loadtxt(file_id, skiprows=3)
            if(spin_occupation==1): #{{{
                print(Header+"** Warning: we attempt to describe splitting of the electrons with different spins.")
                matrixelements_spinned=np.multiply(np.zeros((matrix_side*2, matrix_side*2)), 0e0+0e0j)
                for i in np.arange(0,2*matrix_side):
                    for j in np.arange(0,2*matrix_side):
                        #duplicate each line and each column by blocks of 2. 
                        matrixelements_spinned[i,j]=matrixelements[i/2, j/2]
                print(matrixelements_spinned)
                for i in np.arange(1,2*matrix_side,2):
                    for j in np.arange(1,2*matrix_side,2):
                        #if a line number is even, coefficients d_11 and d_12 are switched (exchange of spins)
                        print(Header+"Terms to be swapped: (i,j)<->(i-1,j): "+str(i)+", "+str(j)+" <-> "+str(i-1)+", "+str(j))
                        Swap(matrixelements_spinned[i,j], matrixelements_spinned[i-1,j]) #swap coeffs of the same line
            #}}}
                        
                        
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
    #print "Info: successfull import of complex-valued dipolar matrix elements for m=-1, m=0, m=1."
    #print "** Warning: we still don't know how to combine these."
    #print "** Warning: therefore one selects only one file among m=-1,0,1 for now!"
    if(verbose==1):
        print("Info: successfull import of real-valued dipolar matrix elements for X direction.")
    return matrixelements, matrix_side
#}}} 
