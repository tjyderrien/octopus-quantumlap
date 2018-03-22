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

## Captures the parameters input from the command line. 
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

## Prepare 0D E(t) bicolor fields for a series of laser parameters. 
def PrepareBiColorFields0D(): 
  ParameterTable = getInput()
  duration1, amplitude1, wavelength1, phase1, polarization1, delay1, duration2, amplitude2, wavelength2, phase2, polarization2, delay2 = ExtractParameters(ParameterTable)

  print "** Preparing E(t) plots..."

  #===================================================
  print "| Plot E(t) for many carrier envelope phases"
  #===================================================

  phase2list = np.linspace(0e0, pi, 6)

  for phase2 in phase2list:
    FigName = str(1E9*float(wavelength1))+"nm-Phase"+str(phase1)+"+"+str(1E9*float(wavelength2))+"nm-Phase"+str(float(phase2)*180e0/pi)+"-Delay"+str(float(delay2)*1E12)+"ps"
    FigTitle = "$\lambda_1=$"+str(1E9*wavelength1)+" nm, $\phi_1=$"+str(phase1*180e0/pi)+", \n $\lambda_2=$"+str(1E9*wavelength2)+" nm, $\phi_2=$"+str(phase2*180e0/pi)+", $\Delta t = $"+str(delay2*1E12)+" ps"

    tmin = -4e0 * duration1
    tmax = 4e0 * duration1

    t = np.linspace(tmin, tmax, 10000)

    Pulse1 = NonChirpedGaussianPulse(t, amplitude1, wavelength1, phase1, duration1, polarization1, delay1)
    Pulse2 = NonChirpedGaussianPulse(t, amplitude2, wavelength2, phase2, duration2, polarization2, delay2)
    Pulse = Pulse1+Pulse2

    PlotPulse(t, Pulse, FigName, FigTitle)

  #==============================================
  print "| Preparing E(t) for many double pulse delays..."
  #==============================================

  phase2=phase1
  delay2list = np.linspace(-2e0*duration1, 2e0*duration2, 6)

  for delay2 in delay2list:
    FigName = str(1E9*wavelength1)+"nm-Phase"+str(phase1)+"+"+str(1E9*wavelength2)+"nm-Phase"+str(phase2*180e0/pi)+"-Delay"+str(delay2*1E12)+"ps"
    FigTitle = "$\lambda_1=$"+str(1E9*wavelength1)+" nm, $\phi_1=$"+str(phase1*180e0/pi)+", \n $\lambda_2=$"+str(1E9*wavelength2)+" nm, $\phi_2=$"+str(phase2*180e0/pi)+", $\Delta t = $"+str(delay2*1E12)+" ps"

    tmin = -4e0 * duration1
    tmax = 4e0 * duration1

    t = np.linspace(tmin, tmax, 10000)

    Pulse1 = NonChirpedGaussianPulse(t, amplitude1, wavelength1, phase1, duration1, polarization1, delay1)
    Pulse2 = NonChirpedGaussianPulse(t, amplitude2, wavelength2, phase2, duration2, polarization2, delay2)
    Pulse = Pulse1+Pulse2

    PlotPulse(t, Pulse, FigName, FigTitle)


## Prepare 2D figures E(t, X) where X is another continuous pulse parameter
# X is sequentially: 
# - wavelength2
# - amplitude2
def PrepareBiColorFields0D_Batch():
  #===========================================================================================
  print "Prepares E(t, wavelength2) figure for given 2-color fields amplitudes, phases, delay"
  #===========================================================================================

  wavelength1 = 800e-9
  duration1 = 30E-15; duration2 = duration1
  phase1 = 0E0; phase2 = 0E0
  amplitude1 = 5E9; amplitude2 = 5E9
  polarization1 = 0E0; polarization2 = polarization1; 
  delay1=0E0; delay2=0E0

  wavelength2 = np.linspace(300e-9, 2000e-9, 100)
  tmin = -1e0 * duration1
  tmax =  1e0 * duration1
  t = np.linspace(tmin, tmax, 10000)

  tt, wavelength22 = np.meshgrid(t, wavelength2)

  print np.shape(tt), np.shape(wavelength22)

  #exit()
  Pulse1 = NonChirpedGaussianPulse(tt, amplitude1, wavelength1, phase1, duration1, polarization1, delay1)
  Pulse2 = NonChirpedGaussianPulse(tt, amplitude2, wavelength22, phase2, duration2, polarization2, delay2)
  Pulse  = Pulse1+Pulse2

  #print wavelength1, amplitude1, phase1
  filename="Bicolor-wavelength"+str(int(wavelength1*1E9))+"-Efield"+str(round(amplitude1*1E-9))+"Vnm-CEP"+str(int((phase1*180./pi)))

  plt.figure()
  plt.contourf(1E15*tt, 1E6*wavelength22, Pulse, cmap=plt.cm.RdBu_r) #plt.cm.Blues
  plt.xlabel(r'Time t (fs)')
  plt.ylabel(r'$\lambda_2$ (um)')
  plt.colorbar()
  plt.title(r"$\lambda_1=$ "+str(int(wavelength1*1E9))+" nm")
  print "Writing file: "+filename
  plt.savefig(filename+".eps")
  plt.savefig(filename+".png")
  #plt.show()

  del wavelength22

  #===============================================================================================
  print "Prepares E(t, amplitude2) figure for given 2-color fields wavelengths, phases"
  #===============================================================================================
  
  LogAmplitude2 = np.linspace(0., 10., 100)
  amplitude2    = np.power(10., LogAmplitude2)
  wavelength2   = wavelength1 / 2. #1030e-9
  
  tt, amplitude22 = np.meshgrid(t, amplitude2)
  
  Pulse1 = NonChirpedGaussianPulse(tt, amplitude1, wavelength1, phase1, duration1, polarization1, delay1)
  Pulse2 = NonChirpedGaussianPulse(tt, amplitude22, wavelength2, phase2, duration2, polarization2, delay2)
  Pulse  = Pulse1+Pulse2

  #print wavelength1, amplitude1, phase1
  filename="Bicolor-wavelength1-"+str(int(wavelength1*1E9))+"-wavelength2-"+str(int(wavelength2*1E9))+"-Efield"+str(round(amplitude1*1E-9))+"Vnm-CEP"+str(int((phase1*180./pi)))+"-ScanAmplitude2"

  plt.figure()
  plt.contourf(1E15*tt, 1E-9*amplitude22, Pulse, cmap=plt.cm.RdBu_r) #plt.cm.Blues
  plt.xlabel(r'Time t (fs)')
  plt.ylabel(r'$E_2$ (V/nm)')
  plt.colorbar()
  plt.title(r"$\lambda_{1}=$ "+str(int(wavelength1*1E9))+" nm, "+r"$\lambda_{2}=$ "+str(int(wavelength2*1E9))+" nm")
  print "Writing file: "+filename
  plt.savefig(filename+".eps")
  plt.savefig(filename+".png")
  #plt.show()

  del amplitude22
  
  #===============================================================================================
  print "Prepares E(t, TimeDelay) figure for given 2-color fields wavelengths, amplitudes, phases"
  #===============================================================================================
  
  amplitude2    = amplitude1
  delay2 = np.linspace(-duration1, duration1, 100)
  wavelength2   = 1030e-9 # wavelength1 / 2. #
  
  tt, delay22 = np.meshgrid(t, delay2)
  
  Pulse1 = NonChirpedGaussianPulse(tt, amplitude1, wavelength1, phase1, duration1, polarization1, delay1)
  Pulse2 = NonChirpedGaussianPulse(tt, amplitude2, wavelength2, phase2, duration2, polarization2, delay22)
  Pulse  = Pulse1+Pulse2

  #print wavelength1, amplitude1, phase1
  filename="Bicolor-wavelength1-"+str(int(wavelength1*1E9))+"-wavelength2-"+str(int(wavelength2*1E9))+"-Efield1-"+str(round(amplitude1*1E-9))+"-Efield-"+str(round(amplitude1*1E-9))+"Vnm-CEP"+str(int((phase1*180./pi)))+"-ScanDelay2"

  plt.figure()
  plt.contourf(1E15*tt, 1E15*delay22, Pulse, cmap=plt.cm.RdBu_r) #plt.cm.Blues
  plt.xlabel(r'Time t (fs)')
  plt.ylabel(r'$\Delta t$ (fs)')
  plt.colorbar()
  plt.title(r"$\lambda_{1}=$ "+str(int(wavelength1*1E9))+" nm, "+r"$\lambda_{2}=$ "+str(int(wavelength2*1E9))+" nm")
  print "Writing file: "+filename
  plt.savefig(filename+".eps")
  plt.savefig(filename+".png")
  #plt.show()

  del delay22
  
  #===============================================================================================
  print "Prepares E(t, CEP2) figure for given 2-color fields wavelengths, amplitudes, durations"
  #===============================================================================================
  print pi 
  #exit()
  amplitude2    = amplitude1
  delay2 = 0.
  wavelength2   = 1030e-9 #wavelength1 / 2.
  phase2 = np.linspace(-pi, pi, 100)
  
  tt, phase22 = np.meshgrid(t, phase2)
  
  Pulse1 = NonChirpedGaussianPulse(tt, amplitude1, wavelength1, phase1, duration1, polarization1, delay1)
  Pulse2 = NonChirpedGaussianPulse(tt, amplitude2, wavelength2, phase22, duration2, polarization2, delay2)
  Pulse  = Pulse1+Pulse2

  #print wavelength1, amplitude1, phase1
  filename="Bicolor-wavelength1-"+str(int(wavelength1*1E9))+"-wavelength2-"+str(int(wavelength2*1E9))+"-Efield1-"+str(round(amplitude1*1E-9))+"-Efield-"+str(round(amplitude1*1E-9))+"Vnm-CEP"+str(int((phase1*180./pi)))+"-ScanPhase2"

  plt.figure()
  plt.contourf(1E15*tt, phase22*180./pi, Pulse, cmap=plt.cm.RdBu_r) #plt.cm.Blues
  plt.xlabel(r'Time t (fs)')
  plt.ylabel(r'$\phi (deg)$')
  plt.colorbar()
  plt.title(r"$\lambda_{1}=$ "+str(int(wavelength1*1E9))+" nm, "+r"$\lambda_{2}=$ "+str(int(wavelength2*1E9))+" nm")
  print "Writing file: "+filename
  plt.savefig(filename+".eps")
  plt.savefig(filename+".png")
  #plt.show()

  del tt, phase22


#def BuildSinglePulseDB(wavelengthMin, wavelengthMax, numpoints, phaseMin, phaseMax, amplitudeMin, amplitudeMax, durationMin, durationMax, polarization): #{{
  ## Prepare a number of plots and input lines for single pulse database
  
##}}

# PrepareBiColorFields0D()
PrepareBiColorFields0D_Batch()