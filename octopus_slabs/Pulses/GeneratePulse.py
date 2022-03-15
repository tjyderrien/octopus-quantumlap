#!python
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2022 T. J.-Y. Derrien
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

## @package GeneratePulse
# This module defines various types of laser pulses and calls the spp_extended_theory/Keldysh/libKeldysh package to construct a pulse envelope for Octopus TD-DFT package.
# A routine provides the pulse envelope in the octopus:td_from_file-compatible format. 
# See octopus-code.org. 

import sys

from scipy import integrate

from octopus_slabs.Libs.libAtomicUnits import *
from spp_extended_theory.Libs.libDatabase import ExportToTxt
from spp_extended_theory.Keldysh.libKeldyshPulses import PulseSmoothTopHatShape
from spp_extended_theory.Libs.libSPP import *

Header="[octopus_slabs: GeneratePulses.py] "
DebugMode = 0 #0: outputs the envelope of the field for octopus. 1: output the fields and plots. 
## Generates the required pulse envelope from spp-ext-theory/Keldysh/libKeldyshPulses and returns it in a file. 
# Input parameters 
# @param t (SI): time instants to prepare the pulse
# @param wavelength (SI)
# @param pulse duration (SI)
# @param peak field amplitude (SI)
# @param CEP (SI)
# @param t0 (SI): central instant of the pulse to generate
# @param PulseDelay (SI): temporal delay for this pulse (useful to generate several pulses)
# Outputs:
# t (At. Units) | Envelope of the field A(t)
def GeneratePulse(t, tau, tauRise, PeakField, wavelength, CEP=0., t0=0., PulseDelay=0.,power=2): #{{{
  # Tested and validated with octopus on Jan 30th 2018. 
  #Envelope, Field = PulseSquaredSinTemporalShape(t, tau, PeakField, wavelength, CEP, t0, PulseDelay)
  # BUG: To be validated. 
  Envelope, Field = PulseSmoothTopHatShape(t, tau, tauRise, PeakField, wavelength, CEP, t0, PulseDelay, power)
  # We obtain SI electric field envelope. We need to convert it to Vector potential envelope in Atomic Units. 
  # NOTE: Are envelope of the electric field and vector potential the same thing? E(t) = - 1/c * dA/dt 
  # NOTE: A=-i c E/omega <=> E = i A * omega / c
  c_SI = c #libUnits
  omega_SI = 2.*pi*c_SI/wavelength
  omega_AU = 1./(Time_SI_to_AU(1./omega_SI))
  
  # Still exact
  print(Header+"** Input: Peak field    E amplitude (SI): "+str(np.abs(Field).max()))
  print(Header+"** Input: Peak envelope E amplitude (SI): "+str(np.abs(Envelope).max()))
  Field_E_AU = Field_SI_to_AU(Field) #complex
  Env_E_AU   = Field_SI_to_AU(Envelope) #Units from V/m to Atomic Units (eV/Bohr?).
  print(Header+"** Converted: Peak field    E amplitude (at.u.): "+str(np.abs(Field_E_AU).max()))
  print(Header+"** Converted: Peak envelope E amplitude (at.u.): "+str(np.abs(Env_E_AU).max()))
  t_AU = Time_SI_to_AU(t) #VALIDATED
  
  Field_A_AU_approx = 1j * c_AU * Field_E_AU / omega_AU
  Env_A_AU_approx = np.sqrt(Field_A_AU_approx * np.conj(Field_A_AU_approx)) #might solve the problem
  
  if(DebugMode == 0): #approximation with error on phase
    # This approximation is inducing a quadrature in the thing. 
    # This kills the resolution at the cycle. 
    print(Header+"** Converted with approximation: Vector field   A amplitude (at.u.): "+str(Field_A_AU_approx.min())+", "+str(Field_A_AU_approx.max()))
    print(Header+"** Converted with approximation: Envelope field A amplitude (at.u.): "+str(np.abs(Env_A_AU_approx).max()))
    
  elif(DebugMode == 3): #attempt to perform the integration directly
  # This should be more exact IF WE WOULD HAVE A NUMERICAL INTEGRATION METHOD FOR COMPLEX INTEGRALS. 
  #PhaseShift=40./180. * pi
    Field_E_AU_real = Field_E_AU.real
    print(Header+"** Real part of E(t): Peak field E amplitude (at.u.): "+str(Field_E_AU_real.min())+", "+str(Field_E_AU_real.max()))
    #Field_A_AU = - c_AU * Field_E_AU_real.cumsum() #leads to a too high amplitude
    Field_A_AU_real = - c_AU * integrate.cumtrapz(Field_E_AU_real, t_AU, initial=0) #should be done in atomic units, and using a real part
    #Env_A_AU = np.sqrt(Field_A_AU * np.conj(Field_A_AU))
    #Field_A_AU = - c_AU * integrate.romb(Field_E_AU) #should be done in atomic units!
    #Field_A_AU = Field_A_AU * np.exp(-1j*PhaseShift) #correcting the phase to obtain a more smooth envelope
    print(Header+"** Converted with integration: Vector field A amplitude (at.u.): "+str(Field_A_AU_real.min())+", "+str(Field_A_AU_real.max()))
    #print Header+"** Converted with integration: Envelope field A amplitude (at.u.): "+str(Field_A_AU_real.max())
    
  # Extraction of the envelope
  
  # Debug mode output
  if(DebugMode==0): 
    Output1 = t_AU
    Output2 = Env_A_AU_approx.real #real-valued using the analytical approximation
    #Output2 = Field_A_AU_approx #DEBUG PURPUSE ONLY
  elif(DebugMode==1): #to solve the initial phase problem and the continuity of the field
    Output1 = t_AU #[:len(t_AU)-1]
    Output2 = Field_A_AU #complex-valued. ONLY FOR DEBUG ! Otherwise, it may send a sinusoidal ENVELOPE to octopus...
  elif(DebugMode==2): #to solve the peak amplitude problem
    Output1 = t_AU
    Output2 = Field_E_AU #AU units, complex-valued (contains the phase)
  elif(DebugMode==3): #We now check the envelope final shape without using the analytical approximation
    Output1 = t_AU #[:len(t_AU)-1]
    #Output2 = Field_A_AU.real #DEBUG PURPUSE ONLY #Env_A_AU.real #real-valued
    Output2 = Field_A_AU_real #real-valued
  else: 
    print(Header+"** ERROR: debug mode variable is undefined. ")
    sys.exit()
  
  return Output1, Output2, 0.*t_AU #GOOD for octopus if Output is an envelope
  #}}}

if __name__ == "__main__": #{{{

  # Read pulse data from command line
  # =======================================================
  if(len(sys.argv)<7):
    print("Usage: ./GeneratePulse.py \ ")
    print("    <Pulse duration (s)> \ ")
    print("    <Ramp cycle number> \ ")
    print("    <Peak field amplitude (V/m)> \ ")
    print("    <Laser wavelength (m)> \ ")
    print("    <CEP [multiples of Pi] > \ ")
    print("    <Double pulse Delay (s) > \ ")
    print("    <End of simulation time (s) > \ ")
    print("    <Power of sin for increase/decrease: 2, 4, 6...> ")
    
    print("Example: ./GeneratePulse.py 10e-15 1e-15 9e9 800e-9 0. 0. 20e-15 2")
    print("Warning: peak amplitude may not be conserved for pulses close to and lower than 10 fs. ")
    exit()

  #============ Manage the command line input =================
  tau        = sys.argv[1]
  RiseCycle  = sys.argv[2]
  PeakField  = sys.argv[3]
  wavelength = sys.argv[4]
  CEP        = sys.argv[5]
  PulseDelay = sys.argv[6]
  ComputingT = sys.argv[7] #maximum time of computation
  power      = sys.argv[8] #sin2 or sin4 for the Soft top hat

  tau = np.asfarray(tau); 
  RiseCycle = np.asfarray(RiseCycle); 
  PeakField = np.asfarray(PeakField);   
  wavelength = np.asfarray(wavelength); 
  CEP = np.asfarray(CEP); 
  PulseDelay = np.asfarray(PulseDelay); 
  ComputingT = np.asfarray(ComputingT);
  power      = np.asfarray(power)     ;

  omega = 2.*pi*c/wavelength

  tauRise = 2.*pi*RiseCycle / (omega)

  tauSafe = wavelength / c #1 optical cycle to ensure having dA/dt = 0 at instant zero !

  ComputingT_min = 7./5.*(tauRise+tau+tauRise+PulseDelay) #this will avoid runs to terminate too early
  if(ComputingT_min > ComputingT):
    print(Header+"** Error: minimum computer time is too small (currently: "+str(ComputingT)+"). Use minimum: "+str(ComputingT_min))
#    print(Header+"** TODO: replace the total value of computing time in inp as well.")
#    exit()
    print(Header+"** Warning: the value was rewritten to avoid computing for nothing.")
    ComputingT = ComputingT_min


  #print Header+"Pulse duration: "+str(tau)

  print(Header+"Computing time = "+str(float(ComputingT)))
  tmin = 0.
  tmax = float(ComputingT)
  Number = 50000 #NOTE: better to put too many nodes here ! Better to have something wider than planned octopus time steps. 
  t = np.arange(tmin, tmax, (tmax-tmin)/(Number+1.))
  #print t
  #exit()
  Envelope     = GeneratePulse(t, tau, tauRise, PeakField, wavelength, 2.*pi*CEP, tau/2.+tauRise+tauSafe, PulseDelay+0.1E-15, power) #NOTE: we always add 0.1 fs to ensure a zero derivative at initial time

  if(DebugMode==1):
    time, Ax, Ay = Envelope
    print(time, Ax)
    plt.figure()
    plt.plot(time, Ax/Ax.max(), label=r"Norm. $A(t)$")
    differential = np.diff(Ax)
    plt.plot(time[:len(time)-1], differential/differential.max(), label=r"Norm. $\partial A / \partial t$") #NOTE: this should be zero in instant 0 to avoid jumps in current! dA/dt should also be continuous. 
    plt.legend()
    plt.grid()
    plt.show()
      
  elif(DebugMode==2):
    time, Ex, Ey = Envelope #AU
    print(time, Ex) #complex-valued
    
    omega_AU = 1./(Time_SI_to_AU(1./omega))
    
    #Ax_cumsum  = -1.j * c_AU * Ex.cumsum() #AU, WARNING. cumsum is not applicable on a complex-valued quantity !!
    Ax_quad    = - c_AU * integrate.cumtrapz(Ex, time) #AU 
    Ax_approx  = c_AU * Ex / omega_AU #NOTE: modest analytical approximation
    #Ax_peak   # has no meaning! 
    
    print(np.shape(time))
    
    plt.figure()
    plt.title(r"$\tau=$"+str(round(tau*1E15))+r" fs, $\tau_{rise}=$"+str(round(tauRise*1E15))+" fs.")
    #print np.shape(Ax_approx)
    plt.plot(time, np.abs(Ax_approx), label=r"$|A(t)|$ approx.")
    plt.plot(time, Ax_approx.real, label=r"Re$(A(t))$")
    #print np.shape(Ax_cumsum)
    #plt.plot(time, Ax_cumsum,  label=r"$A(t)$ cumsum")
    #print np.shape(Ax_quad)
    plt.plot(time[:len(time)-1], np.abs(Ax_quad),  label=r"$|A(t)|$ quad")
    plt.plot(time[:len(time)-1], Ax_quad.real,  label=r"Re$(A(t))$ quad")
    #differential = np.diff(Ex)
    #plt.plot(time[:len(time)-1], differential/differential.max(), label=r"Norm. $\partial A / \partial t$") #NOTE: this should be zero in instant 0 to avoid jumps in current! But should it be continuous? 
    plt.legend()
    plt.grid()
    plt.show()    

  elif(DebugMode==3):
    time, Ax, Ay = Envelope
    print(Header+"** Info: Max value of Ax = "+str(Ax.max()))
    plt.figure()
    plt.plot(time, Ax, label=r"$A(t)$")
    plt.legend()
    plt.grid()
    plt.show()

  OctopusReady = np.transpose(Envelope)

  #print OctopusReady

  ExportToTxt(OctopusReady, "Pulses.dat")

  print(Header+"Successfully generated the pulse. See Pulses.dat.")



#}}}
