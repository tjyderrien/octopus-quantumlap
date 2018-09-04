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


## @module libPlottingSettings
## Define settings for matplotlib
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
from matplotlib.legend_handler import HandlerLine2D

rc('mathtext', default='regular')
rc('text', usetex=False)
rcParams['font.family']='serif'
rcParams['font.size']  = 14
rcParams['font.serif'] =['Times-Roman', 'Bitstream Charter', 'Palatino', 'DejaVu Serif']

# rcParams['font.family']     = 'sans-serif'
# rcParams['font.sans-serif'] = ['phvl8a'] #, 'Helvetica', 'DejaVu Sans']
#['cmtt10'] is failling on negative signs. 

# rc('font',**{'family':'sans-serif','sans-serif':['cm'], 'size':'14'})
# for Palatino and other serif fonts use:
# rc('font', **{'family':'serif', 'serif':['Palatino'], 'size':'14'})
mp.rcParams['legend.numpoints'] = 1

