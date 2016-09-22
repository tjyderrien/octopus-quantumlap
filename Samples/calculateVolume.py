#!/usr/bin/env python
#-*- coding: utf-8 -*-

#IMPORT LIBRARIES
import numpy as np
from numpy import genfromtxt, loadtxt, chararray

from scipy.constants import c, epsilon_0, mu_0, pi, e, m_e, h
import sys

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
  Material = loadtxt(Filename, dtype='str', delimiter=' ', skiprows=2)

  # Calcule the volume of the cell
  X = np.asfarray(Material[:,1])
  Y = np.asfarray(Material[:,2])
  Z = np.asfarray(Material[:,3])

  volume = ( X.max() - X.min() ) * ( Y.max() - Y.min() ) * ( Z.max() - Z.min() )

  return volume #}}}
