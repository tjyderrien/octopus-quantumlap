#!/usr/bin/env python2
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2020 T. J.-Y. Derrien
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

#IMPORT LIBRARIES
import sys
import os
import glob
import itertools

from libDatabase import * #uses symbolic link
#from ../../../spp-extended-theory/SimpleSPP/libDatabase import FilterDatabaseContains

## Returns the name of the most recent octopus logfile
def LookUpLatestLogFile(Filename): #{{{{
  Header = "[libString: LookUpLatestLogFile] "
  import glob
  newest = max(glob.iglob(Filename+'*.log'), key=os.path.getctime)
  print Header+"List of files: "+str(newest)
  # We shall split the path by slashes. 
  newest_splitted = newest.split("/")
  print Header+"List of files: "+str(newest_splitted)
  newest_file = newest_splitted[len(newest_splitted)-1] #last element of the path is selected
  print Header+"List of files: "+str(newest_splitted)
  result = [int(s) for s in newest_file.split(".") if s.isdigit()]
  print Header+"List of files: "+str(result)
  # For some folders, result was empty (?!)
  try:
    final = max(result)
  except: 
    final = result
  return final
#}}}

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
		print(os.path.join(dirname, subdirname))

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
    lists = sorted(itertools.izip(*[names, numbers, directions]))
    names_sorted, numbers_sorted, directions_sorted = list(itertools.izip(*lists))
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
    

## Capture eigenvalues for a given k-point from "info" file
# @param kpoints_index: index of the k-point of interest
# @param filename: "info" file.
def GetEigenValue(kpoints_index, filename="info", verbose=1):
    Header="[libString: GetEigenValue(): ]"
    pattern= "#k ="
    if(verbose==1): 
        print "Info: Construction of the eigen-values matrix."
        print "Pattern of interest: \""+pattern+"\""
        
    with open( filename, "r" ) as source:
        for heading_and_lines in group_by_heading( source, pattern ):
            heading = heading_and_lines[0]
            lines = heading_and_lines[1:] #BUG ? It seems we select 1 k point too far in the list
            if(str(kpoints_index)) in heading:
                #print heading
                #print lines
                break #only keep the first
    print Header+str(heading)
    eigenvalues=[]
    for id in lines:
        line = id.split()
        eigenvalues.append(line[2])
    #print eigenvalues
    return np.asfarray(eigenvalues)

## Among a list of files, extract the one that matches with the Gamma point
# @param filename: common part of the multiple filenames
# @param kpoints: number of kpoints in the simulation
# @listoffiles: full list of the dipolar files
# @CaptureDirectGap: bool for automatic detection.
def ExtractFileForGivenKpoint(filename, kpoints, listoffiles, CaptureDirectGap=True, verbose=False, field_polarization_dir=0): #{{{
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
    expected_num_files = spatial_directions * kpoints**3
    if(len(listoffiles) != expected_num_files):
        print "Error: set of data may be incomplete. Expected number of files: "+str(expected_num_files)
    else: 
        if(verbose==1): 
            print "Info: number of files is correct."
        else: 
            print ""

    # 3. Construct the big matrice of dipolar transition inside the Floquet Hamiltonian. 
    # WARNING: how to put m=-1, m=0, m=+1 ? does this order matter? 

    # 3.1: NOTE: we start with gamma point only. 
    # Gamma point index may be indicated in "info" file, at the gap indication. 
    
    #print "** WARNING: this part is experimental!"
    ## METHOD 1: trusting the choice of octopus
    info_file = "info"
    
    print "Capturing the Gamma point."
    # NOTE: eigen values are sorted by norm. Then taking the first k point should be sufficient. 
    if(CaptureDirectGap):
        print "Is there a direct gap value?"    
        try: 
            kpoints_index_inp   = LookUpStringInFile(info_file, "Direct gap")
            for i in np.arange(0,10):
                kpoints_index_inp=kpoints_index_inp.replace("  ", " ")
                if(verbose==1):
                    print Header+"** Info: decreasing number of spaces"
                    print kpoints_index_inp
            kpoints_index_array = kpoints_index_inp.split(" ")
            #print kpoints_index_array
            kpoints_index       = int(kpoints_index_array[4])
            gamma_point_files_index=[kpoints_index]
            print "** Info: found a direct gap: "+str(gamma_point_files_index)
        except:
            print "** ERROR: could not find a direct gap. Looking for an indirect gap then."
            exit()
            #kpoints_index_inp   = LookUpStringInFile(info_file, "Indirect gap")
            #kpoints_index_array = kpoints_index_inp.split(" ")
            #kpoints_index       = int(kpoints_index_array[7]) #BUG: dangerous...
            #gamma_point_files_index=[kpoints_index]
            ##gamma_point_files_index=[kpoints_index*spatial_directions, kpoints_index*spatial_directions+1, kpoints_index*spatial_directions+2]
    else:
        print "Selecting the first k-point, corresponding to Gamma."
        #kpoints_index_inp   = LookUpStringInFile(info_file, "k #1")
        #for i in np.arange(0,10):
            #kpoints_index_inp=kpoints_index_inp.replace("  ", " ")
            #if(verbose==1):
                #print Header+"** Info: decreasing number of spaces"
                #print kpoints_index_inp
        #kpoints_index_array = kpoints_index_inp.split(" ")
        ##print kpoints_index_array
        print "** Warning: kpoints_index selection is MANUAL."
        kpoints_index       = -1
        gamma_point_files_index=[kpoints_index]
    # Managing filenames for dipolar matrix elements
    gamma_point_files = []
    for num in gamma_point_files_index:
        gamma_point_file = filename+'.k'+str(num)+directions[field_polarization_dir]
        gamma_point_files.append(gamma_point_file)
    if(verbose==1):
        print "Files containing the matrix at Gamma point: "+str(gamma_point_files)
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
        print "Info: size of coupling matrix elements: "+str(matrix_side)+"x"+str(matrix_side)

    for file_id in gamma_point_files: #
        if(verbose==1): 
            print "Info: opening file "+file_id
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
                print Header+"Info: opening "+file_id
            matrixelements=np.loadtxt(file_id, skiprows=3)
            if(spin_occupation==1): #{{{
                print Header+"** Warning: we attempt to describe splitting of the electrons with different spins."
                matrixelements_spinned=np.multiply(np.zeros((matrix_side*2, matrix_side*2)), 0e0+0e0j)
                for i in np.arange(0,2*matrix_side):
                    for j in np.arange(0,2*matrix_side):
                        #duplicate each line and each column by blocks of 2. 
                        matrixelements_spinned[i,j]=matrixelements[i/2, j/2]
                print matrixelements_spinned
                for i in np.arange(1,2*matrix_side,2):
                    for j in np.arange(1,2*matrix_side,2):
                        #if a line number is even, coefficients d_11 and d_12 are switched (exchange of spins)
                        print Header+"Terms to be swapped: (i,j)<->(i-1,j): "+str(i)+", "+str(j)+" <-> "+str(i-1)+", "+str(j)
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
        print "Info: successfull import of real-valued dipolar matrix elements for X direction."
    return matrixelements, matrix_side
#}}} 
