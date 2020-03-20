#!/usr/bin/env python2
#-*- coding: utf-8 -*-

# Copyright (C) 2016-2020 T. J.-Y. Derrien
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

# @module plotPumpProbe
# Provides a set of functions aimed at analyzing the optical properties of solids from the total_current. 

import numpy as np
from libAtomicUnits      import *
from libPlotting         import *
from libPlottingSettings import *

SizeX = 8
SizeY = 12
plotA = 211; plotB=212

ShrinkBy=0

# 1. Select 2 runs,
path_run1_folder="lambda800e-9-tau20E-15-field1E10-CEP0.0-SinSquared2/lambda400e-9-tau20E-15-field0E0-CEP0.0-delay0.5e-15-SinSquared-InducedFields/td.general/"
path_run1 = "lambda800e-9-tau20E-15-field1E10-CEP0.0-SinSquared2/lambda400e-9-tau20E-15-field0E0-CEP0.0-delay0.5e-15-SinSquared-InducedFields/td.general/total_current" #PUMP 
path_run2 = "lambda800e-9-tau20E-15-field1E10-CEP0.0-SinSquared2-Probe1.0-ProbeDelay1e-15/lambda400e-9-tau20E-15-field0E0-CEP0.0-delay0.5e-15-SinSquared-InducedFields/td.general/total_current" #PUMP+PROBE RUN

# 2. Capture their total_current. 
filename = path_run1_folder+"../inp"
current1=np.loadtxt(path_run1, dtype='str', skiprows=4)
current2=np.loadtxt(path_run2, dtype='str', skiprows=4)
# Considered using ASE here, but it is poorly developed for TDDFT. 

# 3. Identify the common temporal intervals. 
current1_time_AU, current1_X_AU, current1_Y_AU, current1_Z_AU = GetTotalCurrent(current1)
current2_time_AU, current2_X_AU, current2_Y_AU, current2_Z_AU = GetTotalCurrent(current2)

current12_time_AU_min = np.max([np.min(current1_time_AU), np.min(current2_time_AU)])
current12_time_AU_max = np.min([np.max(current1_time_AU), np.max(current2_time_AU)])

# 4. Interpolate the signals if necessary (seems not required here). 
# 5. Prepare Yabana's case: J(pump+probe) - J_pump = J_probe. 
current_probe_X_AU=current2_X_AU-current1_X_AU
current_probe_Y_AU=current2_Y_AU-current1_Y_AU
current_probe_Z_AU=current2_Z_AU-current1_Z_AU

plt.figure(figsize=(SizeX,SizeY))
ax1 = plt.subplot(plotA)
plt.ylabel(r'$J_{tot}(t)$ (arb. u.)')
plot_current_PumpProbe, = ax1.plot(current1_time_AU, current2_X_AU, 'g-', label=r'$J_{P+p}(t)$ (a.u.)')
plot_current_Pump,      = ax1.plot(current1_time_AU, current1_X_AU, 'r-', label=r'$J_P(t)$ (a.u.)')
plot_current_probe,     = ax1.plot(current1_time_AU, current_probe_X_AU, 'b-', label=r'$J_p(t)$ (a.u.)')
plot1   = [plot_current_PumpProbe, plot_current_Pump, plot_current_probe]
labels1 = [l.get_label() for l in plot1]

# 6. Fourier Transform the J_probe. 
PulseDuration_SI = GetParameterValueFromInputFile(filename, "tauAsi = ")
wavelength_SI    = GetParameterValueFromInputFile(filename, "lambdaA = ")

current_probe_X_AU_damped = SignalDampR(Time_AU_to_SI(current1_time_AU), current_probe_X_AU,  wavelength_SI, 3*PulseDuration_SI) 

maxSize=current1_time_AU.size-ShrinkBy

TimeSteps=np.diff(current1_time_AU)
TimeStep=TimeSteps[0]

HHG_AU=np.fft.fft(current_probe_X_AU_damped)
HHG_AU_omega=2.*pi*np.fft.fftfreq(maxSize,d=TimeStep*Time_AU_to_SI(1.))
# c_SI = Velocity_AU_to_SI(c_AU)
HHG_AU_order=HHG_AU_omega/(2.*pi*c/wavelength_SI)
HHG_AU_duration=1./HHG_AU_omega

ax2 = plt.subplot(plotB)
Xaxis = HHG_AU_omega[:maxHHG]*hbar/e #energy of the harmonics
Yaxis = np.power(HHG_AU_order[:maxHHG],2.)*np.abs(HHG_AU[:maxHHG])**2
plt.plot(Xaxis, Yaxis, '-', label=r'HHG$(\omega)$') #normalization in amplitude
plt.legend(loc="best")
plt.show()
# 7. Extract the conductivity out of it. 
# 8. Deduce the dielectric property. 
