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

## @package GeneratePulse
# This module defines various types of laser pulses and calls the spp-extended-theory/Keldysh/libKeldysh package to construct a pulse envelope for octopus. 
# A routine provides the pulse envelope in an octopus:td_from_file-compatible format. 

from libKeldyshPulses import PulseSmoothTopHatShape
from libDatabase      import ExportToTxt
from libSPP           import *
from libAtomicUnits   import *

Header="[octopus-slabs: GeneratePulses.py] "

## Generates the required pulse envelope from spp-ext-theory/Keldysh/libKeldyshPulses and returns it in a file. 
# Input parameters 
# @param t (SI): time instants to prepare the pulse
# @param wavelength (SI)
# @param pulse duration (SI)
# @param peak field amplitude (SI)
# @param CEP (SI)
# @param t0 (SI): time zero of the origin
# @param PulseDelay (SI): temporal delay for this pulse (useful to generate several pulses)
# Outputs:
# t (At. Units) | Envelope of the field A(t)
def GeneratePulse(t, tau, tauRise, PeakField, wavelength, CEP=0., t0=0., PulseDelay=0.): #{{{
  # Tested and validated with octopus on Jan 30th 2018. 
  #Envelope, Field = PulseSquaredSinTemporalShape(t, tau, PeakField, wavelength, CEP, t0, PulseDelay)
  # BUG: To be validated. 
  Envelope, Field = PulseSmoothTopHatShape(t, tau, tauRise, PeakField, wavelength, CEP, t0, PulseDelay)
  # We obtain SI electric field envelope. We need to convert it to Vector potential envelope in Atomic Units. 
  # NOTE: Are envelope of the electric field and vector potential the same thing? E(t) = - 1/c * dA/dt 
  # NOTE: A=-i c E/omega <=> E = i A * omega / c
  c_SI = c #libUnits
  omega_SI = 2.*pi*c_SI/wavelength
  omega_AU = 1./(Time_SI_to_AU(1./omega_SI))
  Env_E_AU = Field_SI_to_AU(Envelope) #Units from V/m to Atomic Units (eV/Bohr?). 
  Env_A_AU = c_AU * Env_E_AU / omega_AU
  t_AU = Time_SI_to_AU(t) #VALIDATED
  return t_AU, Env_A_AU, 0.*t_AU
  #}}}

# Read pulse data from command line
# =======================================================
if(len(sys.argv)<6):
  print "Usage: ./GeneratePulse.py \ "
  print "    <Pulse duration (s)> \ "
  print "    <Duration of the ramp (s)> \ "
  print "    <Peak field amplitude (V/m)> \ "
  print "    <Laser wavelength (m)> \ "
  print "    <CEP [multiples of Pi] > \ "
  print "    <Double pulse Delay (s) > \ "
  print "    <End of simulation time (s) > "
  
  print "Example: ./GeneratePulse.py 10e-15 1e-15 9e9 800e-9 0. 0. 20e-15 "
  print "Warning: peak amplitude may not be conserved for pulses close to and lower than 10 fs. "
  exit()

#============ Manage the command line input =================
tau        = sys.argv[1]
tauRise    = sys.argv[2]
PeakField  = sys.argv[3]
wavelength = sys.argv[4]
CEP        = sys.argv[5]
PulseDelay = sys.argv[6]
ComputingT = sys.argv[7] #maximum time of computation

tau = np.asfarray(tau); tauRise = np.asfarray(tauRise); 
PeakField = np.asfarray(PeakField);   wavelength = np.asfarray(wavelength); 
CEP = np.asfarray(CEP); PulseDelay = np.asfarray(PulseDelay); 
ComputingT = np.asfarray(ComputingT); 

#print Header+"Pulse duration: "+str(tau)

print Header+"Computing time = "+str(float(ComputingT))
tmin = 0.
tmax = float(ComputingT)
Number = 200
t = np.arange(tmin, tmax, (tmax-tmin)/(Number+1.))
#print t
#exit()
Envelope     = GeneratePulse(t, tau, tauRise, PeakField, wavelength, CEP, tau, PulseDelay)
OctopusReady = np.transpose(Envelope)

#print OctopusReady

ExportToTxt(OctopusReady, "Pulses.dat")

print Header+"Successfully generated the pulse. See Pulses.dat."



