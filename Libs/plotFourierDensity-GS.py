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
import matplotlib as mp
import matplotlib.pyplot as plt

# from pylab import *

# Settings for matplotlib
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':'16'})
## for Palatino and other serif fonts use:
#rc('font', **{'family':'serif', 'serif':['Palatino'], 'size':'18'})
#rc('text', usetex=True)
mp.rcParams['legend.numpoints'] = 1

from scipy.io import netcdf
filecontent = netcdf.netcdf_file('density.ncdf', 'r')

density = filecontent.variables['rdata'][0,:,:]

# Let's plot these data to know their orientation. 
plt.figure()
plt.xlabel('X')
plt.ylabel('Y')
plt.contourf(density)
plt.show()