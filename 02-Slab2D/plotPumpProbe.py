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

## @package plotPumpProbe
# Provides a set of functions aimed at analyzing the optical properties of solids from the total_current. 

import numpy as np
from libAtomicUnits      import *
from libPlotting         import *
from libPlottingSettings import *

SizeX = 8
SizeY = 12
plotA = 211; plotB=212
LineWidthAll=0.5

ShrinkBy=0
maxHHG = 10*200 #maximum number of samples to cut the spectrum with

# 1. Select 2 runs,

# PUMP PULSE
path_run1 = "lambda800e-9-tau20E-15-field1E10-CEP0.0-SinSquared2/lambda400e-9-tau20E-15-field0E0-CEP0.0-delay0.5e-15-SinSquared-InducedFields/td.general/" 

current1=np.loadtxt(path_run1+"total_current", dtype='str', skiprows=4)


current1_time_AU, current1_X_AU, current1_Y_AU, current1_Z_AU = GetTotalCurrent(current1)
current1_time_SI = Time_AU_to_SI(current1_time_AU)

maxSize=current1_time_AU.size-ShrinkBy

filename = path_run1+"../inp"
PulseDuration_SI = GetParameterValueFromInputFile(filename, "tauAsi = ")
wavelength_SI    = GetParameterValueFromInputFile(filename, "lambdaA = ")

TimeSteps=np.diff(current1_time_AU)
TimeStep=TimeSteps[0]

# Plot the corresponding PUMP current
# (plotting the probe current is too small: useless)
plt.figure(figsize=(SizeX,SizeY))
ax1 = plt.subplot(plotA)
plt.xlabel('Time (s)')
plt.ylabel(r'Pump current $J_P(t)$ (arb. u.)')
plot_current_Pump,      = ax1.plot(current1_time_SI, current1_X_AU, 'r-', label=r'$J_P(t)$ (a.u.)')
# plot1   = [plot_current_PumpProbe, plot_current_Pump, plot_current_probe]
plot1     = [plot_current_Pump]

ax2   = plt.subplot(plotB)
plot2 = []

ax12 = ax1.twinx()
ax12.set_ylabel("Probe current $J_p(t)$ (a.u.)")

# PUMP+PROBE PULSES
for delay in ["1e-15", "5e-15", "10e-15", "15e-15", "30e-15"]:
    #path_run2 = "lambda800e-9-tau20E-15-field1E10-CEP0.0-SinSquared2-Probe1.0-ProbeDelay1e-15/lambda400e-9-tau20E-15-field0E0-CEP0.0-delay0.5e-15-SinSquared-InducedFields/td.general/" #PUMP+PROBE RUN
    path_run2="lambda800e-9-tau20E-15-field1E10-CEP0.0-SinSquared-SinPower2-2-Probe1E-3-ProbeDelay"+str(delay)+"/lambda400e-9-tau20E-15-field0E0-CEP0.0-delay0.5e-15-SinSquared-InducedFields/td.general/" 

    # 2. Capture their total_current. 
    current2=np.loadtxt(path_run2+"total_current", dtype='str', skiprows=4)
    # Considered using ASE here, but it is poorly developed for TDDFT. 

    # 3. Identify the common temporal intervals. 
    current2_time_AU, current2_X_AU, current2_Y_AU, current2_Z_AU = GetTotalCurrent(current2)
    current12_time_AU_min = np.max([np.min(current1_time_AU), np.min(current2_time_AU)])
    current12_time_AU_max = np.min([np.max(current1_time_AU), np.max(current2_time_AU)])

    # 4. Interpolate the signals if necessary (seems not required here). 
    # 5. Prepare Yabana's case: J(pump+probe) - J_pump = J_probe. 
    current_probe_X_AU=current2_X_AU-current1_X_AU
    current_probe_Y_AU=current2_Y_AU-current1_Y_AU
    current_probe_Z_AU=current2_Z_AU-current1_Z_AU

    #plot_current_PumpProbe, = ax2.plot(current1_time_SI, current2_X_AU, '--', label=r'$J_{P+p}(t)$ (a.u.)')
    #plot_current_probe,     = ax2.plot(current1_time_SI, current_probe_X_AU, '-', label=r'$J_p(t)$ (a.u.)')

    # 6. Fourier Transform the J_probe. 
    current_probe_X_AU_damped = SignalDampR(current1_time_SI, current_probe_X_AU,  wavelength_SI, 3*PulseDuration_SI)
    plot_current_probe, = ax12.plot(current1_time_SI, current_probe_X_AU_damped, '--', linewidth=LineWidthAll, label=r'$J_{p}(t),$ $\Delta t =$'+str(np.round(1E15*float(delay),0))+' fs')
    plot1+=[plot_current_probe]
    
    HHG_AU=np.fft.fft(current_probe_X_AU_damped)
    HHG_AU_omega=2.*pi*np.fft.fftfreq(maxSize,d=TimeStep*Time_AU_to_SI(1.))
    # c_SI = Velocity_AU_to_SI(c_AU)
    HHG_AU_order=HHG_AU_omega/(2.*pi*c/wavelength_SI)
    HHG_AU_duration=1./HHG_AU_omega

    Xaxis = HHG_AU_omega[:maxHHG]*hbar/e #energy of the harmonics
    Yaxis = np.power(HHG_AU_order[:maxHHG],2.)*np.abs(HHG_AU[:maxHHG])**2
#    Yaxis_Re = np.power(HHG_AU_order[:maxHHG],2.)*np.real(HHG_AU[:maxHHG])
#    Yaxis_Im = np.power(HHG_AU_order[:maxHHG],2.)*np.imag(HHG_AU[:maxHHG])
    plot_spectrum_probe, = ax2.semilogy(Xaxis, Yaxis, '-', linewidth=LineWidthAll, label=r'$J^2(\omega),$ $\Delta t =$'+str(np.round(1E15*float(delay),0))+' fs') #normalization in amplitude
#    plot_spectrum_probe_Im, = ax2.semilogy(Xaxis, Yaxis_Im, '--', linewidth=LineWidthAll, label=r'$Im(J(\omega)),$ $\Delta t =$'+str(np.round(1E15*float(delay),0))+' fs')
    plot2 += [plot_spectrum_probe]

labels1 = [l.get_label() for l in plot1]
ax1.legend(labels1, loc="best")

labels2 = [l.get_label() for l in plot2]
ax2.set_xlabel("Photon energy (eV)")
ax2.set_ylabel("Harmonic ampl. (arb. u.)")
ax2.legend(labels2, loc="best")
# plt.savefig("Current_Fourier_omega-delay-"+str(float(delay)*1E15)+"fs.eps")
plt.tight_layout()
plt.savefig("Current_Fourier_omega-MultipleDelays.eps")
plt.show()

# 7. Fourier transform of the E_field


# 8. Extract the conductivity out of it. 
# Fourier(J)(omega) = Conductivity(omega) x Fourier ( Field ) ( omega )

# 9. Deduce the dielectric property. 
