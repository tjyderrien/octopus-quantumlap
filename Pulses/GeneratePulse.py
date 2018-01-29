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
# This module defines various types of laser pulses and calls the spp-extended-theory/Keldysh/libKeldysh package. 

from libKeldyshPulses import PulseSquaredSinTemporalShape
from libDatabase      import ExportToTxt
from libSPP           import *

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
# Output
# t (At. Units) | Envelope of the field #NOTE: field or vector potential? Difference? 
def GeneratePulse(t, tau, PeakField, wavelength, CEP, t0, PulseDelay): #{{{
  Envelope, Field = PulseSquaredSinTemporalShape(t, tau, PeakField, wavelength, CEP=0., t0=0., PulseDelay=0.)
  return t, Envelope
  #}}}

# Read pulse data from command line
# =======================================================
if(len(sys.argv)<6):
  print "Usage: ./GeneratePulse.py           \ "
  print "    <Pulse duration> \ "
  print "    <Peak field amplitude (V/m)> \ "
  print "    <Laser wavelength (m)>"
  print "Example: ./GeneratePulse.py 10e-15 9e9 800e-9 0. 0. 20e-15 "
  exit()

#============ Manage the command line input =================
tau        = sys.argv[1]
PeakField  = sys.argv[2]
wavelength = sys.argv[3]
CEP        = sys.argv[4]
#t0         = sys.argv[5]
PulseDelay = sys.argv[5]
ComputingT = sys.argv[6] #maximum time of computation

tau = np.asfarray(tau); PeakField = np.asfarray(PeakField);   wavelength = np.asfarray(wavelength); 
CEP = np.asfarray(CEP); PulseDelay = np.asfarray(PulseDelay); ComputingT = np.asfarray(ComputingT); 

print Header+"Computing time = "+str(float(ComputingT))
tmin = 0.
tmax = float(ComputingT)
t = np.arange(tmin, tmax, (tmax-tmin)/200.)
#print t
#exit()
Envelope     = GeneratePulse(t, tau, PeakField, wavelength, CEP, 2.*tau, PulseDelay=0.)
OctopusReady = np.transpose(Envelope)
ExportToTxt(OctopusReady, "Pulses.dat")



