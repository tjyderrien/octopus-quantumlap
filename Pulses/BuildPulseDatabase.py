#!/usr/bin/env python
#-*- coding: utf-8 -*-

from libPulseGenerator import *

lengthunit=1e0

# =========================================================================
# Set of parameters for generating the 8-dimensions pulse database
# =========================================================================
# Pulse 1
# - Wavelength 1
# - Pulse duration
# - Phase inside the envelope
# - Polarization direction
# - Peak amplitude
# Pulse 2
# - Wavelength 2
# - Pulse duration 2
# - Phase inside the envelope 2 
# - Polarization direction 2
# - Peak amplitude 2
# - Temporal delay with pulse 1. 

ConditionOnInput = (len(sys.argv) != 12)

def getInput(): #{{
  if(ConditionOnInput):
    print "Usage: ./file.py <duration1 amplitude1 wavelength1 phase1 polarization1 delay1> <duration2 amplitude2 wavelength2 phase2 polarization2 delay2>"
    print ""
    print "Using default arguments: "
    
    # For a first example, we take. 
    duration1=10e-15; amplitude1 = 1e9; wavelength1 = 1030e-9; phase1 = 0e0; polarization1="s"; delay1=0e0
    duration2=10e-15; amplitude2 = 1e9; wavelength2 = 515e-9; phase2 = pi/3e0; polarization2="s"; delay2=0e0
    
    ParameterTable=np.zeros((2,6))
    
  ParameterTable=np.array([[duration1, amplitude1, wavelength1, phase1, polarization1, delay1], [duration2, amplitude2, wavelength2, phase2, polarization2, delay2]])
  printCurrentParameters(ParameterTable)
  return ParameterTable
#}}

# TODO: put coefficients in arguments / or in input file, to call the script from bash in scriptable way. 
# TODO: Export the math function required to model the pulse. 
# TODO: In a second stage, we can simply output the field map in complex space, and to introduce it directly to octopus. 

ParameterTable = getInput()
duration1, amplitude1, wavelength1, phase1, polarization1, delay1, duration2, amplitude2, wavelength2, phase2, polarization2, delay2 = ExtractParameters(ParameterTable)

#exit()

print "** PLOTTING..."

#==============================================
print "| Effect of the phase..."
#==============================================

phase2list = np.linspace(0e0, pi, 6)

for phase2 in phase2list:
  FigName = str(1E9*float(wavelength1))+"nm-Phase"+str(phase1)+"+"+str(1E9*float(wavelength2))+"nm-Phase"+str(float(phase2)*180e0/pi)+"-Delay"+str(float(delay2)*1E12)+"ps"
  FigTitle = "$\lambda_1=$"+str(1E9*wavelength1)+" nm, $\phi_1=$"+str(phase1*180e0/pi)+", \n $\lambda_2=$"+str(1E9*wavelength2)+" nm, $\phi_2=$"+str(phase2*180e0/pi)+", $\Delta t = $"+str(delay2*1E12)+" ps"

  tmin = -4e0 * duration1
  tmax = 4e0 * duration1

  t = np.linspace(tmin, tmax, 10000)

  Pulse1 = NonChirpedPulse(t, amplitude1, wavelength1, phase1, duration1, polarization1, delay1)
  Pulse2 = NonChirpedPulse(t, amplitude2, wavelength2, phase2, duration2, polarization2, delay2)
  Pulse = Pulse1+Pulse2

  PlotPulse(t, Pulse, FigName, FigTitle)

#==============================================
print "| Effect of the double pulse delay..."
#==============================================

phase2=phase1
delay2list = np.linspace(-2e0*duration1, 2e0*duration2, 6)

for delay2 in delay2list:
  FigName = str(1E9*wavelength1)+"nm-Phase"+str(phase1)+"+"+str(1E9*wavelength2)+"nm-Phase"+str(phase2*180e0/pi)+"-Delay"+str(delay2*1E12)+"ps"
  FigTitle = "$\lambda_1=$"+str(1E9*wavelength1)+" nm, $\phi_1=$"+str(phase1*180e0/pi)+", \n $\lambda_2=$"+str(1E9*wavelength2)+" nm, $\phi_2=$"+str(phase2*180e0/pi)+", $\Delta t = $"+str(delay2*1E12)+" ps"

  tmin = -4e0 * duration1
  tmax = 4e0 * duration1

  t = np.linspace(tmin, tmax, 10000)

  Pulse1 = NonChirpedPulse(t, amplitude1, wavelength1, phase1, duration1, polarization1, delay1)
  Pulse2 = NonChirpedPulse(t, amplitude2, wavelength2, phase2, duration2, polarization2, delay2)
  Pulse = Pulse1+Pulse2

  PlotPulse(t, Pulse, FigName, FigTitle)

#def BuildSinglePulseDB(wavelengthMin, wavelengthMax, numpoints, phaseMin, phaseMax, amplitudeMin, amplitudeMax, durationMin, durationMax, polarization): #{{
  ## Prepare a number of plots and input lines for single pulse database
  
##}}
