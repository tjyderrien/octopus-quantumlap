#!/usr/bin/env python2
#-*- coding: utf-8 -*-

# Copyright (C) 2016-2019 T. J.-Y. Derrien
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

import numpy as np
from numpy import genfromtxt, loadtxt, chararray, fft, heaviside
# from pylab import *
import sys
import decimal

from libAtomicUnits import *
from libPlottingSettings import *
from libPlotting import *

SizeX = 10
SizeY = SizeX #/GoldenNumber

maxHHG = 3*200 #maximum number of samples to cut the spectrum with

Header="[octopus-slabs: plotLinearResponse] "

TotalCurrent=loadtxt("total_current", dtype='str', skiprows=4)

timeCurrent=np.asfarray(TotalCurrent[:,1]) #Unit a.u.
CurrentX=np.asfarray(TotalCurrent[:,2])

wavelength_SI=1.

# CurrentX_HHG = SignalDampR(Time_AU_to_SI(timeCurrent), CurrentX,  wavelength_SI, tMaxSimulation_SI)

maxSize=301; TimeStep=1

HHG_AU=fft.fft(CurrentX) #NOTE: not damped yet
HHG_AU_omega=2.*pi*fft.fftfreq(maxSize,d=TimeStep*Time_AU_to_SI(1.))
# c_SI = Velocity_AU_to_SI(c_AU)
HHG_AU_order=HHG_AU_omega/(2.*pi*c/wavelength_SI)
HHG_AU_duration=1./HHG_AU_omega

plot_current, = plt.plot(timeCurrent, CurrentX, 'r-', label=r'$J(t)$ (a.u.)')
# plot_current_damped, = plt.plot(timeCurrent, CurrentX_HHG, 'b--', label=r'$J(t)$ (a.u.)')
Xaxis = HHG_AU_omega[:maxHHG]*hbar/e #energy of the harmonics
Yaxis = np.power(HHG_AU_order[:maxHHG],2.)*np.abs(HHG_AU[:maxHHG])**2
plt.plot(Xaxis, Yaxis, '-', label=r'HHG$(\omega)$') #normalization in amplitude
plt.xticks(np.arange(1,30, step=2))
plt.show()
plt.savefig("FFT-LinearResponse.eps")
plt.savefig("FFT-LinearResponse.png")

