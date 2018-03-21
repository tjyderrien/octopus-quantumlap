#!/usr/bin/env python2
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2018 T. J.-Y. Derrien
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

# @module calculateVolume
# Computes the volume of computational cell. Limited to orthogonal cells for now! 

#IMPORT LIBRARIES
import numpy as np
from numpy import genfromtxt, loadtxt, chararray

from scipy.constants import c, epsilon_0, mu_0, pi, e, m_e, h
import sys
import os

from libString import *

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
  newest = max(glob.iglob(Filename+'octopus*.log'), key=os.path.getctime)

  #Open the log file
  searchfile = open(newest, "r") 

  # Search into this file the line with volume/cell
  for line in searchfile:
    if "volume/point" in line: 
      result = line; 
  VolPerPoint = float(result[62:75])
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
    print "** Error: VacuumThickness_inp(): Cannot found PeriodicDimensions in inp."
    exit()
  return VacuumThickness
#}}}

## Computes XY area of the sample
def SampleAreaXY(Filename):#{{{
  Material = loadtxt(Filename+"Si.xyz", dtype='str', delimiter=' ', skiprows=2)

  # Calcule the volume of the cell
  X = np.asfarray(filter(None,Material[:,1]))
  Y = np.asfarray(filter(None,Material[:,2]))
  Z = np.asfarray(filter(None,Material[:,3]))

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

def VolumeOfSupercell_Log(Filename): #{{{
  """      Get the total volume of the super cell from the most recent log file
  """

# Select the latest log file

# Get the local content. 
# dirlist = os.listdir(Folder)
# from pprint import pprint
# pprint(dirlist)

# Find there the latest LOG files.
  import glob
  newest = max(glob.iglob(Filename+'octopus*.log'), key=os.path.getctime)

  #Open the log file
  searchfile = open(newest, "r")

  # Search into this file the line with volume/cell
  for line in searchfile:
    if "Cell volume" in line:
      result = line;
  SuperCellVolume = float(result[15:35])
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
  X = np.asfarray(filter(None, Material[:,1]))
  Y = np.asfarray(filter(None, Material[:,2]))
  Z = np.asfarray(filter(None, Material[:,3]))

  volume = ( X.max() - X.min() ) * ( Y.max() - Y.min() ) * ( Z.max() - Z.min() )

  return volume #}}}
