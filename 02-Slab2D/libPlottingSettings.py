#!/usr/bin/env python
#-*- coding: utf-8 -*-

## @module libPlottingSettings
## Define settings for matplotlib
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.legend_handler import HandlerLine2D

# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':'16'})
## for Palatino and other serif fonts use:
# rc('font', **{'family':'serif', 'serif':['Palatino'], 'size':'18'})
rc('text', usetex=False) #True: Does not work on Draco. Just put False, then, and it will be smooth. 
mp.rcParams['legend.numpoints'] = 1

