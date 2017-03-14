#!/usr/bin/env python
#-*- coding: utf-8 -*-

## @module libPlotting
## Define shared functions for plotting routines

from numpy import genfromtxt, loadtxt, chararray
import numpy as np

## Load geometry, n_el, laser and total_current files (atomic units)
def LoadRunOutputs(path):
  geometry = path+"../Si.xyz" #static method
  Ne=loadtxt(path+"../td.general/n_ex", dtype='str', skiprows=9) #total number of excited electrons, not a density!
  #laser=loadtxt("../td.general/laser", dtype='str', skiprows=6) #should be replaced by genfromtxt method, as some errors are present in laser file
  laser=np.genfromtxt(path+"../td.general/laser", dtype='str', skip_header=9, loose=False, invalid_raise=True, filling_values=0.)
  # PostTreatedDensity = loadtxt("postCalculationNexcited.log", dtype='str', skiprows=1) #Unit: number of electrons, which should be multiplied by cell volume (AU). 
  TotalCurrent=loadtxt(path+"../td.general/total_current", dtype='str', skiprows=4)
  Energy = np.genfromtxt(path+"../td.general/energy", dtype='str', skip_header=5)
  return geometry, Ne, laser, TotalCurrent, Energy

