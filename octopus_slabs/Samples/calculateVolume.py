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

## @package calculateVolume
# Computes the volume of computational cell. Limited to orthogonal cells for now! 

import os

# IMPORT LIBRARIES
import numpy as np
from octopus_slabs.Libs.libString import *
from numpy import loadtxt

Header = "[calculateVolume.py] "

def VolumePerPoint_Log(Filename): #{{{
  """      Get the volume per point from the most recent log file
  """
  
# Select the latest log file

# Get the local content. 
# dirlist = os.listdir(Folder)
# from pprint import pprint
# pprint(dirlist)

# Find there the latest LOG files.
  import glob
  newest = max(glob.iglob(Filename+'octopus*.log'), key=os.path.getmtime)

  #Open the log file
  try:
    searchfile = open(newest, "r", encoding="utf-8", errors='ignore') 
  except:
    logging.warning("Could not open the latest Octopus LOG file.")
  
  # Search into this file the line with volume/cell
  for line in searchfile:
    if 'volume/point' in line: 
      position_start = line.find('volume/point')
      result = line[position_start:];
      #print Header+"** Debug info: "+str(result) #valid
  print(Header+"** Info: filtering the right hand side of the matching string: "+str(result))
  result_filtered = result
  result_filtered = result_filtered.replace('volume/point','')
  result_filtered = result_filtered.replace('=','')
  result_filtered = result_filtered.replace('[b^3]','')
  print(Header+"** Info: filtering the volume/point expression from LOGfile: "+str(result_filtered))
  VolPerPoint = float(result_filtered)
#    else:
#      print "** Error: Volume could not be captured in LOG file."
#      volume = 0.
  return VolPerPoint
#}}}

## Gets the volume of vacuum from octopus input file
def VacuumThickness_inp(Filename): #{{{

  try: 
    result = LookUpStringInFile(Filename+"inp", "Lvac = ")
  except: 
    result = 0.
  IsBulkString = LookUpStringInFile(Filename+"inp", "PeriodicDimensions")
  if(float(IsBulkString[20:23])==2): #Slab case
    VacuumThickness = 2.*float(result[6:20]) #Bohr #1 per side of sample!
  elif(float(IsBulkString[20:23])==3): ##Bulk case
    VacuumThickness = 0.
  else: 
    print("** Error: VacuumThickness_inp(): Cannot found PeriodicDimensions in inp.")
    exit()
  return VacuumThickness
#}}}

## Computes XY area of the sample
def SampleAreaXY(Filename):#{{{
  Material = loadtxt(Filename+"Si.xyz", dtype='str', delimiter=' ', skiprows=2)

  # Calcule the volume of the cell
  X = np.asfarray([_f for _f in Material[:,1] if _f])
  Y = np.asfarray([_f for _f in Material[:,2] if _f])
  Z = np.asfarray([_f for _f in Material[:,3] if _f])

  area = ( X.max() - X.min() ) * ( Y.max() - Y.min() )

  return area
#}}}

## Returns the volume of vacuum area (Bohr**3)
def VolumeOfVacuum(Filename): #{{{
  areaXY=SampleAreaXY(Filename)
  thickness=VacuumThickness_inp(Filename)
  return areaXY*thickness
#}}}

  
#}}}

## Get the total volume of the super cell from the most recent log file
# 
def VolumeOfSupercell_Log(Filename): #{{{
  # Select the latest log file

  # Get the local content. 
  # dirlist = os.listdir(Folder)
  # from pprint import pprint
  # pprint(dirlist)

  # Find there the latest LOG files.
  import glob
  newest = max(glob.iglob(Filename+'octopus*.log'), key=os.path.getmtime)
  print(Header+"** Info: Opening LOG file:"+str(newest))
  #Open the log file
  searchfile = open(newest, "r")
  
  # BUG : what if the run just failed? Then a file is defined, but without the relevant information. 

  # Search into this file the line with volume/cell
  for line in searchfile:
    if "Cell volume" in line:
      result = line;
      break;
    else: 
      result = '0'
  print(Header+"** Info: Filtering the volume capture in LOG file.")
  result_filtered = result.replace('Cell volume', '')
  result_filtered = result_filtered.replace('[b^3]','')
  result_filtered = result_filtered.replace(' ', '')
  result_filtered = result_filtered.replace('=', '')
  print(Header+"** Info: Filtered the volume capture.")
  SuperCellVolume = float(result_filtered) #BUG: fails on Cell volume =         90106.3068 [b^3]
#    else:
#      print "** Error: Volume could not be captured in LOG file."
#      volume = 0.
  return SuperCellVolume
#}}}


def VolumeXYZ(Filename): #{{{
  """ Calculate the volume from a XYZ file <Filename>
      Units can be anything. 
  """
  # Loading the XYZ file located into the folder. 
  # Get the filename given in argument. 
  # try: 
  #   Filename = sys.argv[1]
  # except:
  #   print "Usage: ./calculateVolume.py <XYZ file>"
  # exit()

  # Loading the XYZ file located into the folder.
  Material = loadtxt(Filename+"Si.xyz", dtype='str', delimiter=' ', skiprows=2)

  # Calcule the volume of the cell
  X = np.asfarray([_f for _f in Material[:,1] if _f])
  Y = np.asfarray([_f for _f in Material[:,2] if _f])
  Z = np.asfarray([_f for _f in Material[:,3] if _f])

  volume = ( X.max() - X.min() ) * ( Y.max() - Y.min() ) * ( Z.max() - Z.min() )

  return volume #}}}
