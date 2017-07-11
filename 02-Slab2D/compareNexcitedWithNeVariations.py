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

# IMPORT LIBRARIES
import numpy as np
# from numpy import genfromtxt, loadtxt, chararray
# import matplotlib as mp
# import matplotlib.pyplot as plt

# from matplotlib import rc
# from pylab import *
from scipy.constants import c, epsilon_0, mu_0, pi, e, m_e, h
from scipy.integrate import tplquad
# from matplotlib.legend_handler import HandlerLine2D
import sys

# Settings for matplotlib
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':'16'})
## for Palatino and other serif fonts use:
#rc('font', **{'family':'serif', 'serif':['Palatino'], 'size':'18'})
#rc('text', usetex=True)
# mp.rcParams['legend.numpoints'] = 1

from scipy.io import netcdf

# print "Info: getting the GS density."
NeGSdata = netcdf.netcdf_file('../td.0000000/density.ncdf', 'r')
NeGS = NeGSdata.variables['rdata'][:,:,:]
#NeGSdata.close()

# print "Info: getting the TD data."
NeTDdata = netcdf.netcdf_file('density.ncdf', 'r')
NeTD = NeTDdata.variables['rdata'][:,:,:]
#NeTDdata.close()

# print "Info: Building the excess of density..."
#TODO: Use maybe another method: f+ = (|f| - f)/2
NeTDvar = NeTD-NeGS
condlist = [NeTDvar>0.]
choicelist = [NeTDvar]
NeTDpos = np.select(condlist, choicelist)
#NeTDpos2 = 0.5 * ( abs(NeTDvar) + NeTDvar ) #Valid: equal to NeTDpos

# print "Info: Building the depletion of density..."
condlist = [NeTDvar<0.]
choicelist = [NeTDvar]
NeTDneg = np.select(condlist, choicelist)

# Square of the TD-GS to avoid first order increase
NeTD2=NeTDvar**2


# print "Calculate the quantity of excess electrons..."
# Direct sum
TotalNeGS = np.sum(NeGS) #TODO: is volume element missing?
TotalNeTD = np.sum(NeTD) #TODO: is volume element missing? 
TotalNeTDpos = np.sum(NeTDpos) #TODO: is volume element missing?
TotalNeTDneg = np.sum(NeTDneg) #TODO: is volume element missing?
TotalNeTD2 = np.sum(NeTD2) 
# TPLquad integration method
# scipy.integrate.tplquad(NeTDpos)


# print "Post-calculated N_excited:"+str(TotalNeTDpos)
print TotalNeTDpos, TotalNeTDneg, TotalNeTD, TotalNeGS, TotalNeTD2

exit()
