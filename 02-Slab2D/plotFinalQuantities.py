#!/usr/bin/env python
#-*- coding: utf-8 -*-

# @module plotFinalQuantities
# Build the generated octopus simulations
# Builds a table with laser parameters, final densities, and final energies after simulation. 

## IMPORT LIBRARIES
import numpy as np
from numpy import genfromtxt, loadtxt, chararray, fft
# from pylab import *
import sys
import decimal

import calculateVolume as Sample
import libKeldysh as Keldysh
from libAtomicUnits import *
from libPlottingSettings import *
from libString import *

