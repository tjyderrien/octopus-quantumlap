#!/usr/bin/env python
#-*- coding: utf-8 -*-

# IMPORT LIBRARIES
import numpy as np
from numpy import genfromtxt, loadtxt, chararray
import matplotlib as mp
import matplotlib.pyplot as plt

from matplotlib import rc
# from pylab import *
from scipy.constants import c, epsilon_0, mu_0, pi, e, m_e, h
from matplotlib.legend_handler import HandlerLine2D
import sys

# Settings for matplotlib
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':'16'})
## for Palatino and other serif fonts use:
#rc('font', **{'family':'serif', 'serif':['Palatino'], 'size':'18'})
#rc('text', usetex=True)
mp.rcParams['legend.numpoints'] = 1

from scipy.io import netcdf

# print "Info: getting the GS density."
NeGSdata = netcdf.netcdf_file('../td.0000000/density.ncdf', 'r')
NeGS = NeGSdata.variables['rdata'][:,:,:]

# print "Info: getting the TD data."
NeTDdata = netcdf.netcdf_file('density.ncdf', 'r')
NeTD = NeTDdata.variables['rdata'][:,:,:]

# print "Info: Building the excess of density..."
NeTDvar = NeTD-NeGS
condlist = [NeTDvar>0.]
choicelist = [NeTDvar]
NeTDpos = np.select(condlist, choicelist)

# print "Info: Building the depletion of density..."
condlist = [NeTDvar<0.]
choicelist = [NeTDvar]
NeTDneg = np.select(condlist, choicelist)

# print "Calculate the quantity of excess electrons..."
TotalNeTDpos = np.sum(NeTDpos)
# print "Post-calculated N_excited:"+str(TotalNeTDpos)
print TotalNeTDpos

exit()

print "Comparing with the number of excited electrons."
