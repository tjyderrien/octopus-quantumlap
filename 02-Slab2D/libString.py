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

## Splitting, sorting and re-packing the table in a list of filenames.
def GetOrderedListOfFiles(Filename): #{{{
    filelist = glob.iglob(Filename) #makes selection over files, but not sorted
    splitted_t = [[],[]]
    numbers=[]; names=[]
    for file in filelist: #{{{
        splited=file.split(".")
        names.append(splited[0])
        numbers.append(int(splited[1]))
    #}}}
    numbers.sort()
    filelist = np.transpose([names, numbers])

    finallist = []
    for line in filelist:
        newline=str(line[0])+"."+str(line[1])
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
def GetEigenValue(kpoints_index, filename="info"):
    print "Info: Construction of the eigen-values matrix."
    pattern = "#k ="
    print "Pattern of interest: \""+pattern+"\""
    with open( filename, "r" ) as source:
        for heading_and_lines in group_by_heading( source, pattern ):
            heading = heading_and_lines[0]
            lines = heading_and_lines[1:]
            if(str(kpoints_index)) in heading:
                #print heading
                #print lines
                break #only keep the first
    print heading
    eigenvalues=[]
    for id in lines:
        line = id.split()
        eigenvalues.append(line[2])
    #print eigenvalues
    return np.asfarray(eigenvalues)
