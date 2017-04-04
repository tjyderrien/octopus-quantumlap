#!/usr/bin/env python
#-*- coding: utf-8 -*-

## @module Look for strings in files.

#IMPORT LIBRARIES
import sys
import os

from libDatabase import * #uses symbolic link
#from ../../../spp-extended-theory/SimpleSPP/libDatabase import FilterDatabaseContains

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
  
