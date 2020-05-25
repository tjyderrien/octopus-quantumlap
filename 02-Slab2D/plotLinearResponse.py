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

## @package plotLinearResponse
# Provides a set of functions aimed at analyzing the linear optical properties of solids from the total_current. 

import numpy as np
from libAtomicUnits      import *
from libPlotting         import *
from libPlottingSettings import *

SizeX = 12
SizeY = 8
plotA = 121; plotB=122
LineWidthAll=1

EnergyMin=0. #eV
EnergyMax=25 #eV

ShrinkBy=0
maxHHG = 10*200 #maximum number of samples to cut the spectrum with
StopAt = 0.9 #damping starts at 90% of the available timeline
FilterBeforeKick=False

# Switch in plotB
plotEomega=False #< 1E-10 a.u.
plotJomega=False #<1E-2 a.u.
plotConductivity=False #should be max(J)/max(E). Should be 1E3 in experiments.
plotPermittivity=True #should be close to 1E10 in a.u. Strange. Should not depend on unit system.

# 1. Select 2 ruOutput = geometry + density + elf + current + kinetic_energy_density + potential #GENERATES LARGE NUMBER OF FILES.ns,


# PUMP PULSE
# Mo linear response: path_run1 = "lambda800e-9-tau20E-15-field0E0-CEP0.0-SinSquared-SinPower2-2-Probe1E-3-ProbeDelay1e-15/lambda400e-9-tau20E-15-field0E0-CEP0.0-delay0.5e-15-SinSquared-InducedFields/td.general/"
# path_run1   = "lambda800e-9-tau20E-15-field0E0-CEP0.0-SinSquared-SinPower2-2-Probe1E-3-ProbeDelay1e-15/lambda400e-9-tau20E-15-field0E0-CEP0.0-delay0.5e-15-SinSquared-InducedFields/td.general/"
path_run1   = "../td.general/"
current1    =np.loadtxt(path_run1+"total_current", dtype='str', skiprows=4)
gauge_field1=np.loadtxt(path_run1+"gauge_field",   dtype='str', skiprows=4)
# laser       =np.loadtxt(path_run1+"laser",         dtype='str', skiprows=4)

A1tot_time_AU, A1tot_X_AU, A1tot_Y_AU, A1tot_Z_AU             = GetTotalGaugeField(gauge_field1) #BUG? Do we want the probe field here? It is actually equal to A_tot(t)=A_probe*Dirac(t-t0) + A_induced(t)
current1_time_AU, current1_X_AU, current1_Y_AU, current1_Z_AU = GetTotalCurrent(current1)
# Asource_time_AU, Asource_X_AU, Asource_Y_AU, Asource_Z_AU     = #Have to reconstruct the kick signal and to remove it

A1tot_time_SI    = Time_AU_to_SI(A1tot_time_AU)
current1_time_SI = Time_AU_to_SI(current1_time_AU)

maxSize=current1_time_AU.size-ShrinkBy
tMaxSimulation_AU = current1_time_AU[maxSize-1] #this is simulation that has already run
tMaxSimulation_SI = Time_AU_to_SI(tMaxSimulation_AU)

filename = path_run1+"../inp"
PulseDuration_SI = GetParameterValueFromInputFile(filename, "tauAsi = ")
wavelength_SI    = GetParameterValueFromInputFile(filename, "lambdaA = ") #NOTE: meaningless since there is no pump. We iwll rather use the frequency of the probe beam, i.e. a wide distribution.
#omega            = 2.*pi*c/wavelength_SI

TimeSteps            = np.diff(current1_time_AU)
TimeSteps_E_tot_AU   = np.diff(A1tot_time_AU)

TimeStep=TimeSteps[0] #AU
TimeStep_E_tot_AU   = TimeSteps_E_tot_AU[0]

plot1 =[]
plot12=[]
# Plot the corresponding probe current
# (plotting the probe current is too small: useless)
plt.figure(figsize=(SizeX,SizeY))
ax1 = plt.subplot(plotA)
ax1.set_xlabel('Time (s)')
ax1.set_ylabel(r'Probe+induced field $A_p(t)+A_i(t)$ (arb. u.)')
plot_field_probe, = ax1.plot(A1tot_time_SI, A1tot_X_AU, 'b--', linewidth=LineWidthAll, label=r'$A_x^{p}(t)+A_x^{ind}(t)$ (a.u.)')
plot1 += [plot_field_probe] #plot_current_PumpProbe]

ax12 = ax1.twinx()
ax12.set_ylabel("Probed current $J_x^{p}(t)$ (a.u.)")
plot_current_probe, = ax12.plot(current1_time_SI, current1_X_AU, 'r--', linewidth=LineWidthAll, label=r'$J_x^{p}(t)$ (a.u.)')
plot12+=[plot_current_probe]


ax2   = plt.subplot(plotB)
plot2 = []

# PUMP+PROBE PULSES
# for delay in ["1e-15", "5e-15", "10e-15", "15e-15", "30e-15"]: #TODO: generate this list from database
    #path_run2 = "lambda800e-9-tau20E-15-field1E10-CEP0.0-SinSquared2-Probe1.0-ProbeDelay1e-15/lambda400e-9-tau20E-15-field0E0-CEP0.0-delay0.5e-15-SinSquared-InducedFields/td.general/" #PUMP+PROBE RUN
    
current_probe_X_AU=current1_X_AU
current_probe_Y_AU=current1_Y_AU
current_probe_Z_AU=current1_Z_AU

    #plot_current_PumpProbe, = ax2.plot(current1_time_SI, current2_X_AU, '--', label=r'$J_{P+p}(t)$ (a.u.)')
    #plot_current_probe,     = ax2.plot(current1_time_SI, current_probe_X_AU, '-', label=r'$J_p(t)$ (a.u.)')

# 6. Fourier Transform the J_probe. 
current_probe_X_AU_damped = SignalDampR(current1_time_SI, current_probe_X_AU,  wavelength_SI, StopAt*current1_time_SI[len(current1_time_SI)-1]) #tMaxSimulation_SI) #WARNING: we might create a problem here by using SI! #TODO: Change maximum duration to the 3/4th of the way to end of simulation
plot_current_probe_damped, = ax12.plot(current1_time_SI, current_probe_X_AU_damped, 'r-', linewidth=LineWidthAll, label=r'$J_{p}^{(damped)}(t)$')
plot12 += [plot_current_probe_damped]
    
current_probe_Fourier_AU       = np.fft.fft(current_probe_X_AU_damped)
current_probe_Fourier_SI       = Current_AU_to_SI(current_probe_Fourier_AU)
print "Max of J_p (a.u.)="
print np.max(current_probe_Fourier_AU)
print "Max of J_p (S.I.)="
print np.max(current_probe_Fourier_SI)
current_probe_Fourier_AU_omega = 2.*pi*np.fft.fftfreq(maxSize,d=TimeStep) #AU
current_probe_Fourier_SI_omega = 2.*pi*np.fft.fftfreq(maxSize,d=TimeStep*Time_AU_to_SI(1.)) #SI
# c_SI = Velocity_AU_to_SI(c_AU)
# current_probe_Fourier_AU_order=current_probe_Fourier_AU_omega/current_probe_Fourier_AU_omega #meaningless
# current_probe_Fourier_AU_duration=1./current_probe_Fourier_AU_omega

Xaxis = current_probe_Fourier_SI_omega[:maxHHG]*hbar/e #energy of the harmonics (SI here, OK)
# Yaxis = np.power(current_probe_Fourier_AU_omega[:maxHHG],2.)*np.abs(current_probe_Fourier_AU[:maxHHG])
Yaxis = np.abs(current_probe_Fourier_AU[:maxHHG])
#    Yaxis_Re = np.power(current_probe_Fourier_AU_order[:maxHHG],2.)*np.real(current_probe_Fourier_AU[:maxHHG])
# Yaxis = np.real(current_probe_Fourier_AU[:maxHHG])
#    Yaxis_Im = np.power(current_probe_Fourier_AU_order[:maxHHG],2.)*np.imag(current_probe_Fourier_AU[:maxHHG])
if(plotJomega):
  plot_spectrum_probe, = ax2.plot(Xaxis, Yaxis, 'r-', linewidth=LineWidthAll, label=r'$|J(\omega)|$') #normalization in amplitude #OK
#    plot_spectrum_probe_Im, = ax2.semilogy(Xaxis, Yaxis_Im, '--', linewidth=LineWidthAll, label=r'$Im(J(\omega)),$ $\Delta t =$'+str(np.round(1E15*float(delay),0))+' fs')
  plot2 += [plot_spectrum_probe]

# 6. Extracting A_probe from A_tot=A_probe + A_source

GaugeFieldDelay_inp  = LookUpStringInFile(filename, "GaugeFieldDelaySI =")
GaugeFieldDelay_tmp    = GaugeFieldDelay_inp.split("=")
try:
  GaugeFieldDelay      = GaugeFieldDelay_tmp[1] #rhs
except: 
  GaugeFieldDelay=0e0
  print(Header+"Warning: Selection of r.h.s. failed in INP. No data in:"+str(GaugeFieldDelay_tmp))

if(FilterBeforeKick):    
	KickInstant=int(Time_SI_to_AU(float(GaugeFieldDelay))/TimeStep+1) #This removes the signal before the kick
else:
	KickInstant=0 #keeps the kick: no filtering
print(Header+"** Info: Kick was performed at instant #"+str(KickInstant))

A1probe_X_AU    = A1tot_X_AU[KickInstant:] #-Asource_X_AU[:len(A1tot_X_AU)] + A1tot_X_AU
A1probe_Y_AU    = A1tot_Y_AU[KickInstant:] #-Asource_Y_AU[:len(A1tot_Y_AU)] + A1tot_Y_AU
A1probe_Z_AU    = A1tot_Z_AU[KickInstant:] #-Asource_Z_AU[:len(A1tot_Z_AU)] + A1tot_Z_AU

# 7. Converting A_field into E_field
E1probe_time_AU = A1tot_time_AU[KickInstant:] #we keep the same time grid for E and A, but remove everything before the kick starts
E1probe_X_AU    = -1./c_AU * np.divide( np.diff(A1probe_X_AU), np.diff(E1probe_time_AU) ) 
E1probe_Y_AU    = -1./c_AU * np.divide( np.diff(A1probe_Y_AU), np.diff(E1probe_time_AU) )
E1probe_Z_AU    = -1./c_AU * np.divide( np.diff(A1probe_Z_AU), np.diff(E1probe_time_AU) )

E1probe_time_SI = Time_AU_to_SI(1.)*E1probe_time_AU
# Keeping the beginning of the signal
E_probe_X_AU_damped = SignalDampR(E1probe_time_SI[:len(E1probe_time_AU)-1], E1probe_X_AU[:len(E1probe_time_AU)-1], wavelength_SI, StopAt*E1probe_time_SI[len(E1probe_time_SI)-1]) # StopAt*tMaxSimulation_SI) #this does not work well. 
plot_field_probe_damped, = ax1.plot(E1probe_time_SI[:len(E1probe_time_AU)-1], c_AU*E_probe_X_AU_damped, 'b-', linewidth=LineWidthAll, label=r'$c \times E_{p}^{(damped)}(t)$')
plot1  += [plot_field_probe_damped]
plot12 += [plot_current_probe_damped]

# 7. Fourier transform of the E_field
E_probe_Fourier_AU=np.fft.fft(E_probe_X_AU_damped)
E_probe_Fourier_SI=Field_AU_to_SI(E_probe_Fourier_AU)
print "Max of E_p (a.u.)="
print np.max(E_probe_Fourier_AU)
print "Max of E_p (S.I.)="
print np.max(E_probe_Fourier_SI)
E_probe_Fourier_AU_omega=2.*pi*np.fft.fftfreq(maxSize,d=TimeStep)
E_probe_Fourier_SI_omega=2.*pi*np.fft.fftfreq(maxSize,d=TimeStep*Time_AU_to_SI(1.))
# E_probe_Fourier_AU_duration=1./E_probe_Fourier_AU_omega

# 8. Extract the conductivity tensor elements. Let's see if one need symmetries. 
# WARNING: this is based on Ohm law, that must be approximative !!!!
# Fourier(J)(omega) = Conductivity(omega) x Fourier ( Field ) ( omega )
conductivity_AU = np.divide(current_probe_Fourier_AU[:len(E_probe_Fourier_AU)], E_probe_Fourier_AU)
conductivity_SI = np.divide(current_probe_Fourier_SI[:len(E_probe_Fourier_SI)], E_probe_Fourier_SI)

print "Max of conductivity (AU)="
print np.max(conductivity_AU)
print "Max of conductivity (SI)="
print np.max(conductivity_SI)


print "Range of frequencies (AU)="
print E_probe_Fourier_AU_omega
print "Range of frequencies (SI)="
print E_probe_Fourier_SI_omega
# 9. Deduce the dielectric property from conductivity. #returns something close to 1E10 a.u. Convert to SI before!
# permittivity_AU = 1. + 4.*np.pi * 1.j * np.divide(conductivity_AU, current_probe_Fourier_AU_omega[:len(conductivity_AU)]) #atomic units. TODO: this should be done in SI, if I understand right. 
permittivity_AU = 1. + 4.*np.pi * 1.j * np.divide(conductivity_AU, E_probe_Fourier_AU_omega[:len(conductivity_AU)])
permittivity_SI = 1. + 1.j / epsilon_0 * np.divide(conductivity_SI, E_probe_Fourier_SI_omega[:len(conductivity_SI)])

print "Dielectric permittivity (AU)="
print np.abs(np.real(permittivity_AU))
print np.abs(np.imag(permittivity_AU))
print "Dielectric permittivity (SI)="
print np.abs(np.real(permittivity_SI))
print np.abs(np.imag(permittivity_SI))

# 10. Plotting
Xaxis =  current_probe_Fourier_SI_omega[:maxHHG]*hbar/e #energy of the harmonics
X0axis = E_probe_Fourier_SI_omega[:maxHHG]*hbar/e #energy of the harmonics (same result as before?)
Y0axis = np.abs(E_probe_Fourier_AU[:maxHHG])
Yaxis =  np.real(conductivity_AU[:maxHHG])
Y2axis = np.imag(conductivity_AU[:maxHHG])
Y3axis = np.real(permittivity_AU[:maxHHG]) #imag(permittivity) is taken as example in Otobe & Yabana PRB 2012 (Ohm law based calculation). 
Y4axis = np.imag(permittivity_AU[:maxHHG])

if(plotEomega): 
	plot_spectrum_field_probe, = ax2.plot(X0axis, Y0axis, 'k-', linewidth=LineWidthAll, label=r'$|E(\omega)|$')
	plot2 += [plot_spectrum_field_probe]

if(plotConductivity):
	plot_spectrum_conductivity_probe,      = ax2.plot(Xaxis, Yaxis, 'b-', linewidth=LineWidthAll, label=r'$Re(\sigma_{xx})(\omega),$') #normalization in amplitude
	plot_spectrum_conductivity_probe_imag, = ax2.plot(Xaxis, Y2axis, 'r-', linewidth=LineWidthAll, label=r'$Im(\sigma_{xx})(\omega),$') #normalization in amplitude
	plot2 += [plot_spectrum_conductivity_probe, plot_spectrum_conductivity_probe_imag]

if(plotPermittivity):
	plot_spectrum_permittivity_probe, = ax2.plot(Xaxis, Y3axis, 'b--', linewidth=LineWidthAll, label=r'$Re(\varepsilon_{xx})(\omega),$') #normalization in amplitude
	plot_spectrum_permittivity_probe_imag, = ax2.plot(Xaxis, Y4axis, 'r--', linewidth=LineWidthAll, label=r'$Im(\varepsilon_{xx})(\omega),$') #normalization in amplitude
	plot2 += [plot_spectrum_permittivity_probe, plot_spectrum_permittivity_probe_imag]

labels1=   [l.get_label() for l in plot1]
labels12 = [l.get_label() for l in plot12]
labels2 =  [l.get_label() for l in plot2]

ax1.legend( labels1,  loc="upper left" )
ax12.legend(labels12, loc="upper right")

ax1.set_yscale( 'symlog', linthreshy=1E-3)
ax12.set_yscale('symlog', linthreshy=1E-3)
# ax2.set_yscale( 'symlog', linthreshy=1e-12)
ax2.set_xlim(EnergyMin, EnergyMax)

ax2.set_xlabel("Photon energy (eV)")
ax2.set_ylabel("Dielectric permittivity")
ax2.legend(labels2, loc="best")
# plt.savefig("Current_Fourier_omega-delay-"+str(float(delay)*1E15)+"fs.eps")
plt.tight_layout()
plt.savefig("LinearResponse_Current_Fourier_omega.eps")
plt.show()
