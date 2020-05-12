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
plotA = 121; plotB=122
LineWidthAll=0.5

ShrinkBy=0
maxHHG = 10*200 #maximum number of samples to cut the spectrum with

plotEomega=False
plotJomega=False
plotConductivity=True
# 1. Select 2 runs,

# PUMP PULSE
# Mo: 
path_run1 = "lambda800e-9-tau20E-15-field1E10-CEP0.0-SinSquared2/lambda400e-9-tau20E-15-field0E0-CEP0.0-delay0.5e-15-SinSquared-InducedFields/td.general/" 
# Au
# Au: path_run1 = "lambda800e-9-tau20E-15-field1E10-CEP0.0-SinSquared-SinPower2-2/lambda400e-9-tau20E-15-field0E0-CEP0.0-delay0.5e-15-SinSquared-InducedFields/td.general/"
# path_run2 = "lambda800e-9-tau20E-15-field0E0-CEP0.0-SinSquared-SinPower2-2-Probe1E-3-ProbeDelay1e-15/lambda400e-9-tau20E-15-field0E0-CEP0.0-delay0.5e-15-SinSquared-InducedFields/td.general/"

current1    =np.loadtxt(path_run1+"total_current", dtype='str', skiprows=4)
gauge_field1=np.loadtxt(path_run1+"gauge_field", dtype='str', skiprows=4)
laser_field1=np.loadtxt(path_run1+"laser", dtype='str', skiprows=4)

A1tot_time_AU, A1tot_X_AU, A1tot_Y_AU, A1tot_Z_AU             = GetTotalGaugeField(gauge_field1)
A1laser_time_AU, A1laser_X_AU, A1laser_Y_AU, A1laser_Z_AU     = GetTotalGaugeField(laser_field1)
current1_time_AU, current1_X_AU, current1_Y_AU, current1_Z_AU = GetTotalCurrent(current1)

A1tot_time_SI    = Time_AU_to_SI(A1tot_time_AU)
A1laser_time_SI  = Time_AU_to_SI(A1laser_time_AU)
current1_time_SI = Time_AU_to_SI(current1_time_AU)

maxSize=current1_time_AU.size-ShrinkBy

filename = path_run1+"../inp"
PulseDuration_SI = GetParameterValueFromInputFile(filename, "tauAsi = ")
wavelength_SI    = GetParameterValueFromInputFile(filename, "lambdaA = ")
omega            = 2.*pi*c/wavelength_SI

TimeSteps            = np.diff(current1_time_AU)
TimeSteps_E_tot_AU   = np.diff(A1tot_time_AU)
TimeSteps_E_laser_AU = np.diff(A1laser_time_AU)

TimeStep=TimeSteps[0]
TimeStep_E_tot_AU   = TimeSteps_E_tot_AU[0]
TimeStep_E_laser_AU = TimeSteps_E_laser_AU[0]
if(TimeStep_E_tot_AU != TimeStep_E_laser_AU):
    print "** WARNING: timeline of gauge field is different from timeline of laser field"
    exit()
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
for delay in ["1e-15", "5e-15", "10e-15", "15e-15", "30e-15"]: #TODO: generate this list from database
    # Mo: 
    path_run2 = "lambda800e-9-tau20E-15-field1E10-CEP0.0-SinSquared-SinPower2-2-Probe1E-3-ProbeDelay"+str(delay)+"/lambda400e-9-tau20E-15-field0E0-CEP0.0-delay0.5e-15-SinSquared-InducedFields/td.general/" #PUMP+PROBE RUN for Mo
    # Au: path_run2="lambda800e-9-tau20E-15-field1E10-CEP0.0-SinSquared-SinPower2-2-Probe1E-3-ProbeDelay"+str(delay)+"/lambda400e-9-tau20E-15-field0E0-CEP0.0-delay0.5e-15-SinSquared-InducedFields/td.general/" 

    # 2. Capture their total_current. 
    current2=np.loadtxt(path_run2+"total_current", dtype='str', skiprows=4)
    # Considered using ASE here, but it is poorly developed for TDDFT. 

    # 3. Identify the common temporal intervals. 
    current2_time_AU, current2_X_AU, current2_Y_AU, current2_Z_AU = GetTotalCurrent(current2)
    current12_time_AU_min = np.max([np.min(current1_time_AU), np.min(current2_time_AU)])
    current12_time_AU_max = np.min([np.max(current1_time_AU), np.max(current2_time_AU)])

    # 4. TODO: CONVERT THE CURRENT FROM AU TO SI. 
    
    # 5. Prepare Yabana's case: J(pump+probe) - J_pump = J_probe. 
    current_probe_X_AU=current2_X_AU-current1_X_AU
    current_probe_Y_AU=current2_Y_AU-current1_Y_AU
    current_probe_Z_AU=current2_Z_AU-current1_Z_AU
    
    # E(pump+probe+induced) - E(pump) = E(probe+induced)
    # WARNING: is time identical for A(gauge_field) and A(laser)? added a test above
    A1probe_time_AU = A1tot_time_AU
    A1probe_X_AU    = A1tot_X_AU - A1laser_X_AU
    A1probe_Y_AU    = A1tot_Y_AU - A1laser_Y_AU
    A1probe_Z_AU    = A1tot_Z_AU - A1laser_Z_AU
    
    #plot_current_PumpProbe, = ax2.plot(current1_time_SI, current2_X_AU, '--', label=r'$J_{P+p}(t)$ (a.u.)')
    #plot_current_probe,     = ax2.plot(current1_time_SI, current_probe_X_AU, '-', label=r'$J_p(t)$ (a.u.)')

    # 6. Fourier Transform the J_probe. 
    current_probe_X_AU_damped = SignalDampR(current1_time_SI, current_probe_X_AU,  wavelength_SI, 3*PulseDuration_SI)
    plot_current_probe, = ax12.plot(current1_time_SI, current_probe_X_AU_damped, '--', linewidth=LineWidthAll, label=r'$J_{p}(t),$ $\Delta t =$'+str(np.round(1E15*float(delay),0))+' fs')
    plot1+=[plot_current_probe]
    
    Current_Fourier_AU=np.fft.fft(current_probe_X_AU_damped)
    Current_omega_AU=2.*pi*np.fft.fftfreq(maxSize,d=TimeStep*Time_AU_to_SI(1.))
    # c_SI = Velocity_AU_to_SI(c_AU)
    HHG_AU_order=Current_omega_AU/(omega)
    #HHG_AU_duration=1./Current_omega_AU
    if(plotJomega):
      Xaxis = Current_omega_AU[:maxHHG]*hbar/e #energy of the harmonics
      Yaxis = np.power(HHG_AU_order[:maxHHG],2.)*np.abs(Current_Fourier_AU[:maxHHG])**2
#    Yaxis_Re = np.power(HHG_AU_order[:maxHHG],2.)*np.real(Current_Fourier_AU[:maxHHG])
#    Yaxis_Im = np.power(HHG_AU_order[:maxHHG],2.)*np.imag(Current_Fourier_AU[:maxHHG])
      plot_spectrum_Jprobe, = ax2.semilogy(Xaxis, Yaxis, '-', linewidth=LineWidthAll, label=r'$J^2(\omega),$ $\Delta t =$'+str(np.round(1E15*float(delay),0))+' fs') #normalization in amplitude
#    plot_spectrum_Jprobe_Im, = ax2.semilogy(Xaxis, Yaxis_Im, '--', linewidth=LineWidthAll, label=r'$Im(J(\omega)),$ $\Delta t =$'+str(np.round(1E15*float(delay),0))+' fs')
      plot2 += [plot_spectrum_Jprobe]
    
    # 7. Converting A_probe_field into E_field
    # E_tot=-1/c * dA_tot/dt #We keep the induced fields
    #A1tot_time_AU, A1tot_X_AU, A1tot_Y_AU, A1tot_Z_AU = GetTotalGaugeField(gauge_field)
    A1probe_time_SI = Time_AU_to_SI(A1probe_time_AU)
    E1probe_time_SI = A1probe_time_SI #we keep the same time grid
    E1probe_X_AU    = -1./c * np.divide( np.diff(A1probe_X_AU), np.diff(E1probe_time_SI) )
    E1probe_Y_AU    = -1./c * np.divide( np.diff(A1probe_Y_AU), np.diff(E1probe_time_SI) )
    E1probe_Z_AU    = -1./c * np.divide( np.diff(A1probe_Z_AU), np.diff(E1probe_time_SI) )
    # 7. TODO: CONVERT THE FIELD FROM AU TO SI? 
    E_probe_X_AU_damped = SignalDampR(E1probe_time_SI[:len(E1probe_time_SI)-1], E1probe_X_AU[:len(E1probe_time_SI)-1], wavelength_SI, 3*PulseDuration_SI)
    # NOTE: One has to use the PROBE field here.

    # 7. Fourier transform of the E_field
    E_Fourier_AU=np.fft.fft(E_probe_X_AU_damped)
    E_omega_AU  =2.*pi*np.fft.fftfreq(maxSize,d=TimeStep_E_tot_AU*Time_AU_to_SI(1.))
    E_order_AU=E_omega_AU/(omega) #hbar assumed equal to 1 here (atomic units)
    #HHG_AU_duration=1./E_omega_AU
    if(plotEomega):
      Xaxis = E_omega_AU[:maxHHG]*hbar/e #energy of the harmonics
      Yaxis = np.power(E_order_AU[:maxHHG],2.)*np.abs(E_Fourier_AU[:maxHHG])**2
      plot_spectrum_Eprobe, = ax2.semilogy(Xaxis, Yaxis, '--', linewidth=LineWidthAll, label=r'$E^2(\omega),$ $\Delta t =$'+str(np.round(1E15*float(delay),0))+' fs') #normalization in amplitude
      plot2 += [plot_spectrum_Eprobe]
    
    # 8. Extract the conductivity out of it. 
    # Fourier(J_probe)(omega) = Conductivity(omega) x Fourier ( Field_probe ) ( omega )
    # For now, we keep it simple: 
    conductivity = np.divide(Current_Fourier_AU[:maxHHG], E_Fourier_AU[:maxHHG])

    ## 9. Deduce the dielectric property from conductivity. 
#    permittivity = 1. + 1j * conductivity / epsilon_0 / omega #NOTE: is this fine in atomic units?
    permittivity = 1. + 4.*np.pi * 1.j * conductivity / omega #atomic units, [Otobe, Yabana et al., Phys. Rev. B (2016)]
    permittivity_omega_AU  =2.*pi*np.fft.fftfreq(maxSize,d=TimeStep*Time_AU_to_SI(1.))
    permittivity_order_AU  = permittivity_omega_AU / omega #hbar is also missing here
    if(plotConductivity):
      Xaxis = permittivity_omega_AU[:maxHHG]*hbar/e #energy of the harmonics
      Yaxis = np.imag(permittivity[:maxHHG]) #imag(permittivity) is taken as example in Otobe & Yabana PRB 2016. 
      plot_spectrum_permittivity_probe, = ax2.plot(Xaxis, Yaxis, '-.', linewidth=LineWidthAll, label=r'$\varepsilon_{xx}(\omega),$ $\Delta t =$'+str(np.round(1E15*float(delay),0))+' fs') #normalization in amplitude
      plot2 += [plot_spectrum_permittivity_probe]
    
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


