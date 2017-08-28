#!/usr/bin/env python
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2017 T. J.-Y. Derrien
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


## @module Look for strings in files.

#IMPORT LIBRARIES
import sys
import os

from libDatabase import * #uses symbolic link
#from ../../../spp-extended-theory/SimpleSPP/libDatabase import FilterDatabaseContains

## Returns the name of the most recent octopus logfile
def LookUpLatestLogFile(Filename): #{{{{
  import glob
  newest = max(glob.iglob(Filename+'*.log'), key=os.path.getctime)
  result = [int(s) for s in newest.split(".") if s.isdigit()]
  # For some folders, result was empty (?!)
  try:
    final = max(result)
  except: 
    final = result
  return final
#}}}

## Looks for a string in a file. 
# If duplicate files are found, the most recent file is selected.
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