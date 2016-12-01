#!/usr/bin/env python
#-*- coding: utf-8 -*-

## @module Look for strings in files.

#IMPORT LIBRARIES
import sys
import os

## Looks for a string in a file. 
# If duplicate files are found, the most recent file is selected.
def LookUpStringInFile(Filename, string): #{{{
  import glob
  newest = max(glob.iglob(Filename), key=os.path.getctime)

  #Open the log file
  searchfile = open(newest, "r")

  # Search into this file the <line> with <string>
  for line in searchfile:
    if string in line:
      result = line;
  return result
#}}}
