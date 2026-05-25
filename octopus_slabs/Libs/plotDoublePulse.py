#!/usr/bin/env python
#-*- coding: utf-8 -*-

# Copyright (C) 2016-2023 T. J.-Y. Derrien
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

## @module plotDoublePulse
# Plots for two-colour / double-pulse TDDFT runs.
# Functions: plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple,
#            plotEffectOfDoublePulse_Color1_Color2_SingleFieldStrength,
#            plotAdjustedCollisionFrequencies_wavelength

from octopus_slabs.Libs._plotCommon import *  # noqa: F401,F403

# Plots (Ne_final | Wpi | Energy | E**2) (Efield1, Efield2) for pulses of fixed wavelengths, same pulse duration, same CEP. 
# THIS FIGURE IS THE ONE THAT CAN BE COMPARED WITH VLADIMIR ZHUKOV. 
# From this function, we perform a  wish to GENERATE a 1D section in this figure to analyze the role of Stark effect in bicolor absorption.
def plotEffectOfDoublePulse_Efield1_Efield2_SpecificColorCouple(PlottingTable, wavelength1=800e-9, wavelength2=1600e-9, PulseDuration=20e-15, PulseType=3, RiseCycles=2, WhatToPlot=0, Egap_SI=2.56*e, RemoveTunneling=False): #{{{
  LogWpi=False                                                                                                                                                                                                                
  LogPlotXY=False
  InducedFields=0.0 #we dont want to enable this for now
  gammaLimit = 1E-10
  meff = 0.2226 #effective mass of electron [Waroquiers, D. et al, Physical Review B, 2013, 87, 075121]
  CEP2_request=-pi/2.
  print("====== Plot: DOUBLE PULSE, SAME POLARIZATION, SAME PULSE DURATIONS, Ne(E1, E2) =====")
  #Filters: double pulses, same wavelengths, same polarization
  PlottingDoublePulses = PlottingTable #get the full database
  PlottingDoublePulses = libDatabase.FilterDatabaseGreaterThan(PlottingDoublePulses, 1E0, libPlotting.IndexTable['Efield1']) #removes fields where Epeak < 1 V/m. #works
  PlottingDoublePulses = libDatabase.FilterDatabaseGreaterThan(PlottingDoublePulses, 1E0, libPlotting.IndexTable['Efield2']) #removes fields where Epeak2 < 1 V/m to select double pulses.
  PlottingDoublePulses = libDatabase.FilterDatabase(PlottingDoublePulses, float(0.), libPlotting.IndexTable['CrossPolarization']) #removes fields where Polarization1 != Polarization2
  PlottingDoublePulses = libDatabase.FilterDatabase(PlottingDoublePulses, PulseDuration, libPlotting.IndexTable['pulses2']) # Removing fields were pulseDuration1 != pulseDuration2
  PlottingDoublePulses = libDatabase.FilterDatabase(PlottingDoublePulses, wavelength1, libPlotting.IndexTable['wavelengths1']) # Selecting color1
  PlottingDoublePulses = libDatabase.FilterDatabase(PlottingDoublePulses, wavelength2, libPlotting.IndexTable['wavelengths2']) # Selecting color2
  PlottingDoublePulses = libDatabase.FilterDatabase(PlottingDoublePulses, PulseType, libPlotting.IndexTable['PulseType']) # Selecting pulse type
  PlottingDoublePulses = libDatabase.FilterDatabase(PlottingDoublePulses, InducedFields, libPlotting.IndexTable['FieldInduced']) #Removing the induced fields. Yeah because they should not be mixed!
  
  if(CEP2_request >= 0):
      PlottingDoublePulses = libDatabase.FilterDatabaseLowerThan(PlottingDoublePulses,   (CEP2_request)*(1.1), libPlotting.IndexTable['CEP2'])
      PlottingDoublePulses = libDatabase.FilterDatabaseGreaterThan(PlottingDoublePulses, (CEP2_request)*(0.9), libPlotting.IndexTable['CEP2']) #keeps only CEP = 0.
  else:
      PlottingDoublePulses = libDatabase.FilterDatabaseGreaterThan(PlottingDoublePulses, (CEP2_request)*(1.1), libPlotting.IndexTable['CEP2'])
      PlottingDoublePulses = libDatabase.FilterDatabaseLowerThan(PlottingDoublePulses,   (CEP2_request)*(0.9), libPlotting.IndexTable['CEP2']) #keeps only CEP = 0.

  # TODO: we shall clean CEP2 as well, as it must influence the results. 
  
  #PlottingDoublePulses = libDatabase.FilterDatabase(PlottingDoublePulses, CEP2, IndexTable['CEP2']) # Must have same CEP
  
  if(RiseCycles != 0): #0 selects all the possible types
    PlottingDoublePulses = libDatabase.FilterDatabase(PlottingDoublePulses, RiseCycles, libPlotting.IndexTable['RiseTime']) # Selecting Rise Time
  #### 2D PLOTTING Density(Efield1, Efield2) for one given set of laser parameters. For each combination, we produce one plot like this. 
  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq, wavelength2_uniq, CEP2_uniq, Efield2_uniq, pulses2_uniq= libPlotting.CountTheNumberOfLaserParameters(PlottingDoublePulses)
  logger.info(Header+"Total number of datasets: "+str(len(set(Efield1_uniq))*len(set(Efield2_uniq))))

  
  TotalField2D      = np.empty((0,3)) #preparing the 2D plot
  DensityMax2D      = np.empty((0,3)) #preparing the 2D plot
  DensityMax2D_grad = np.empty((0,3)) #preparing the 2D plot
  CurrentMax2D      = np.empty((0,3)) #preparing the 2D plot
  EnergyMax2D       = np.empty((0,3)) #preparing the 2D plot
  WpiMax2D          = np.empty((0,3)) #preparing the 2D plot
  
  for Efield1_i in Efield1_uniq: #{{{ #selecting amplitude1
    PlottingDoublePulses_SelectedField1 = libDatabase.FilterDatabase(PlottingDoublePulses, Efield1_i, libPlotting.IndexTable['Efield1'])
    for Efield2_i in Efield2_uniq: #{{{ #selecting amplitude2
      # Selecting Field2
      PlottingDoublePulses_SelectedField2 = libDatabase.FilterDatabase(PlottingDoublePulses_SelectedField1, Efield2_i, libPlotting.IndexTable['Efield2'])
      
      # Constructing the density table
      runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType_tmp, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingDoublePulses_SelectedField2)
# Modules
      
      #print(Header+"Field1: "+str(Efield1)+" m, Field2: "+str(Efield2)+" m.")
      
      if(len(Ne_final) >= 1): 
        Field1_gamma = libKeldysh.gammaKeldysh(Egap_SI, meff, Efield1[0], wavelengths1[0])
        Field2_gamma = libKeldysh.gammaKeldysh(Egap_SI, meff, Efield2[0], wavelengths2[0])
      
        # Tests if regime is tunneling
        if(RemoveTunneling): 
            Field1_Multiphoton = (Field1_gamma > gammaLimit)
            Field2_Multiphoton = (Field2_gamma > gammaLimit)
        else: #we dont check
            Field1_Multiphoton = True
            Field2_Multiphoton = True
      
        if(Field1_Multiphoton and Field2_Multiphoton): #We filter single-pulses in tunneling regime
            ToBeAddedField  =np.array([Efield1[0], Efield2[0], np.power(Efield1[0]+Efield2[0],2)]) #builds the total peak squared field
            ToBeAdded       =np.array([Efield1[0], Efield2[0], Ne_final[0]])                       #BUG: removes some results, not clean enough!
            ToBeAddedGrad   =np.array([.5*c*epsilon_0*np.power(Efield1[0],2), .5*c*epsilon_0*np.power(Efield2[0],2), Ne_final[0]]) #meant to find d(Nexc)/dI_1,2
            ToBeAddedCurrent=np.array([Efield1[0], Efield2[0], TotalCurrent_final[0]]) #BUG: removes some results, not clean enough!
            ToBeAddedEnergy =np.array([Efield1[0], Efield2[0], ElectronEnergy_max[0]]) #BUG: removes some results, not clean enough!
            if(LogWpi):
                ToBeAddedWpi    =np.array([Efield1[0], Efield2[0], np.log10(W_pi[0])]) #BUG: removes some results, not clean enough!
            else: 
                ToBeAddedWpi    =np.array([Efield1[0], Efield2[0], W_pi[0]]) #BUG: removes some results, not clean enough!
            #print ToBeAdded
            DensityMax2D_grad=np.vstack((DensityMax2D_grad, ToBeAddedGrad)) 
            TotalField2D     =np.vstack((TotalField2D     , ToBeAddedField))
            DensityMax2D     =np.vstack((DensityMax2D     , ToBeAdded))
            CurrentMax2D     =np.vstack((CurrentMax2D     , ToBeAddedCurrent))
            EnergyMax2D      =np.vstack((EnergyMax2D      , ToBeAddedEnergy))
            WpiMax2D         =np.vstack((WpiMax2D         , ToBeAddedWpi))
            #print "Cleaning..."
            del ToBeAdded, ToBeAddedCurrent, ToBeAddedEnergy, ToBeAddedWpi, ToBeAddedField
    #}}}
  #}}}
  #print(Header+" ** Debug: PulseType: "+str(PulseType))
  PulseTypeName = int(np.asfarray(PulseType))
  # Stack is built. Export, show and plot on a heat map. 
  logger.info("** Info: exporting the sub-plot2D: Ne_final(Field1, Field2)...")
  logger.info("** Info: useful to visualize available data with GNUplot. ")
  filename_field   = "Bicolor-"+str(PulseTypeName)+"-Color1-"+str(round(wavelength1*1E9))+"-Color2-"+str(round(wavelength2*1E9))+"-"+str(int(PulseDuration*1e15))+"fs-EPeakWithE1E2"
  filename         = "Bicolor-"+str(PulseTypeName)+"-Color1-"+str(round(wavelength1*1E9))+"-Color2-"+str(round(wavelength2*1E9))+"-"+str(int(PulseDuration*1e15))+"fs-NeFinalWithE1E2"
  filename_gradNexc= "Bicolor-"+str(PulseTypeName)+"-Color1-"+str(round(wavelength1*1E9))+"-Color2-"+str(round(wavelength2*1E9))+"-"+str(int(PulseDuration*1e15))+"fs-NexcWithI1I2"
  filename_current = "Bicolor-"+str(PulseTypeName)+"-Color1-"+str(round(wavelength1*1E9))+"-Color2-"+str(round(wavelength2*1E9))+"-"+str(int(PulseDuration*1e15))+"fs-CurrentFinalWithE1E2"
  filename_energy  = "Bicolor-"+str(PulseTypeName)+"-Color1-"+str(round(wavelength1*1E9))+"-Color2-"+str(round(wavelength2*1E9))+"-"+str(int(PulseDuration*1e15))+"fs-EnergyFinalWithE1E2"
  filename_wpi     = "Bicolor-"+str(PulseTypeName)+"-Color1-"+str(round(wavelength1*1E9))+"-Color2-"+str(round(wavelength2*1E9))+"-"+str(int(PulseDuration*1e15))+"fs-WpiWithE1E2"
  # Sorting the data before storing and plotting
  TotalField2Dsorted = TotalField2D[np.lexsort((TotalField2D[:,2], TotalField2D[:,1], TotalField2D[:,0]))]
  DensityMax2Dsorted = DensityMax2D[np.lexsort((DensityMax2D[:,2], DensityMax2D[:,1], DensityMax2D[:,0]))]
  DensityMax2DgradSorted = DensityMax2D_grad[np.lexsort((DensityMax2D_grad[:,2], DensityMax2D_grad[:,1], DensityMax2D_grad[:,0]))]
  CurrentMax2Dsorted = CurrentMax2D[np.lexsort((CurrentMax2D[:,2], CurrentMax2D[:,1], CurrentMax2D[:,0]))]
  EnergyMax2Dsorted  = EnergyMax2D[np.lexsort((EnergyMax2D[:,2], EnergyMax2D[:,1], EnergyMax2D[:,0]))]
  WpiMax2Dsorted     = WpiMax2D[np.lexsort((WpiMax2D[:,2], WpiMax2D[:,1], WpiMax2D[:,0]))]
  
  del TotalField2D, DensityMax2D, CurrentMax2D, EnergyMax2D, WpiMax2D, DensityMax2D_grad
  
  ## Lazy to build up a new function for that for now.
  unit=1E-9
  PlotGradient=False
  if(WhatToPlot == BicolorE1E2plotDictionnary['Nexc']): 
      QuantityTitle=r'$N_{exc}$ (m$^{-3}$), $(\lambda_1, \lambda_2) = $('+str(round(wavelength1*1E9,0))+', '+str(round(wavelength2*1E9,0))+') nm'
      libDatabase.ExportToTxt(DensityMax2Dsorted, filename+".csv")
  elif(WhatToPlot == BicolorE1E2plotDictionnary['NexcGrad']): #We'll plot grad(Nexc)_Intensity to find number of required photons
      filename=filename_gradNexc
      DensityMax2Dsorted = DensityMax2DgradSorted
      QuantityTitle=r'$N_{exc}$, $(\lambda_1, \lambda_2) = $('+str(round(wavelength1*1E9,0))+', '+str(round(wavelength2*1E9,0))+') nm'
      libDatabase.ExportToTxt(DensityMax2DgradSorted, filename+".csv")
      PlotGradient=True
      unit=1E0 #As we compute gradient, let's keep SI units here. 
  elif(WhatToPlot == BicolorE1E2plotDictionnary['E2']): 
      DensityMax2Dsorted = TotalField2Dsorted; 
      filename = filename_field
      QuantityTitle=r"$E_{tot}^2$ (V/nm)$^2$, $(\lambda_1, \lambda_2) = $("+str(round(wavelength1*1E9,0))+", "+str(round(wavelength2*1E9,0))+") nm."
      libDatabase.ExportToTxt(TotalField2Dsorted, filename+".csv")
  elif(WhatToPlot == BicolorE1E2plotDictionnary['J']): 
      #PlotGradient=False
      DensityMax2Dsorted = CurrentMax2Dsorted; 
      filename = filename_current
      QuantityTitle=r"$J_{end}$ (at.u.), $(\lambda_1, \lambda_2) = $("+str(round(wavelength1*1E9,0))+", "+str(round(wavelength2*1E9,0))+") nm."
      libDatabase.ExportToTxt(CurrentMax2Dsorted, filename+".csv")
  elif(WhatToPlot == BicolorE1E2plotDictionnary['Energy']): 
      #PlotGradient=False
      DensityMax2Dsorted = EnergyMax2Dsorted
      filename = filename_energy
      QuantityTitle=r"$\xi_{end}$ (eV), $(\lambda_1, \lambda_2) = $("+str(round(wavelength1*1E9,0))+", "+str(round(wavelength2*1E9,0))+") nm."
      libDatabase.ExportToTxt(EnergyMax2Dsorted, filename+".csv")
  elif(WhatToPlot == BicolorE1E2plotDictionnary['Wpi']): 
      PlotGradient=True
      DensityMax2Dsorted = WpiMax2Dsorted
      filename = filename_wpi
      QuantityTitle=r'$w_{PI}$ (W$^{-1}$.m$^{-1}$), $(\lambda_1, \lambda_2) = $('+str(round(wavelength1*1E9,0))+', '+str(round(wavelength2*1E9,0))+') nm'
      libDatabase.ExportToTxt(WpiMax2Dsorted, filename+".csv")
  else: 
      PlotGradient=False
      logger.error(Header+"** plot_Zhukov_BicolorMaps: wrong request!")
  
  
  # Shaping the tables
  BrutX = unit*DensityMax2Dsorted[:,0] #Efield1
  BrutY = unit*DensityMax2Dsorted[:,1] #Efield2
  BrutZ = DensityMax2Dsorted[:,2]      #Quantity to plot
  
  ### One need to extract interpolation of BrutZ on a regular mapping. 
  if(not len(BrutZ)==0): 
      RegX, RegY, RegZ = libPlotting.plot2dHeatMap(BrutX, BrutY, BrutZ, QuantityTitle, filename, ShowFinalPlot, PlotGradient)
  else: 
      logger.error((Header+"** Error: Selected database is empty."))
      exit()
  
  ## SPECIAL CASE: preparing figure to support Nadya's idea that BICOLOR induces higher absorption than separed single-colors. 
  ## AbsorbedEnergy(Field or intensity)
  ## Intensity1 
  ## Intensity2
  ## Intensity_1+2 ( t )  = 0.5*c*eps0*|E1(t)+E2(t)|**2 #where is optical index? 
  ## We shall select the peak intensity. 
  ## 
  
  if(WhatToPlot == BicolorE1E2plotDictionnary['Energy']):
    # Back to SI unit, but on a regular mesh. 
    Field1   = RegX/unit #back to V/m. 
    Field2   = RegY/unit #back to V/m. 
    
    # 
    
    ## WARNING: induced fields should be necessary for this computation. 
    # Selecting one single line to plot in these data. 
    plt.figure()
    plt.plot()
  
  #if(WhatToPlot == BicolorE1E2plotDictionnary['Nexc']):
    #print "== Performing a section in the two-dimensional Nexc profile."
    ## We will use interp2d in the Brut data
    
    #DB_Field1 = BrutX
    #DB_Field2 = BrutY
    #DB_Wpi    = BrutZ
    
    #print Header+"Data well imported from DB."
    #print DB_Field1.shape, DB_Field2.shape, DB_Wpi.shape
    ### Attempting to use lists of points for interp2d
    
    #WPI_func = interp2d(DB_Field1, DB_Field1, DB_Wpi, fill_value=-1) #interfaces w_PI with E1, E2 values, as a function in
    
    ##EfieldSweep = np.linspace(1E9,1E10, 10)
    #print WPI_func(1E8, 1E8) ##Weirdly returns something out of bounds
    #exit()
    ## Do we need a 2D matrix interpolation? 
    ## This works for 2D matrix interpolation
    #print Header+"** Preparing interpolation of w_pi(E1,E2)..."
    #xdim = int(np.sqrt(len(DB_Field1))) #NOTE: assumes data are square. True for 800x1600, CEP=0 at least.
    #print xdim
    #Wpi_2D_X  = DB_Field1.reshape(xdim, xdim)
    ##print Wpi_2D_X
    #Wpi_2D_Y  = DB_Field2.reshape(xdim, xdim)
    ##print Wpi_2D_Y
    #Wpi_2D    = DB_Wpi.reshape((xdim, xdim))
    ##print Wpi_2D
    #Wpi_X     = Wpi_2D_X[:,0]
    ##print Wpi_X
    #Wpi_Y     = Wpi_2D_Y[0,:]
    ##print Wpi_Y
    #Wpi_2D_t = np.transpose(Wpi_2D)
    ##print Header+"** Interpolating the w_PI..."
    ##InterpolationOrder=1
    ## WARNING: the 800x1600 files have #1: E2**1 !! and #2: E1**2 ! Therefore we had to transpose the matrix. 
    #WPI_func = interp2d(Wpi_X, Wpi_Y, Wpi_2D_t) #interfaces w_PI with E1, E2 values, as a function in RxR. 
    ##WPI_func = interp2d(Wpi_X, Wpi_Y, Wpi_2D)
    #print Header+"Trying to 2D-interpolate at one (E1,E2) point: "
    #print WPI_func(0., 0.) #**2, FieldEnvelopeNormalized2_CGS**2)
  
  del DensityMax2Dsorted
  
  #print "=== CONSTRUCTING the Zhukov maps from the bicolor Keldysh he provided."
  
#}}}

##### 2D PLOTTING Density(wavelength1, wavelength2) for one single field value? 
#Mode 2: fixed fields, same pulse durations, then plot Ne_final(Color1, Color2)
def plotEffectOfDoublePulse_Color1_Color2_SingleFieldStrength(PlottingTable, FieldStrength, WhatToPlot, PulseDuration=20e-15, Tolerance=0.01, Egap_SI=2.56*e): #{{{
  print("====== Plot: DOUBLE PULSE, SAME POLARIZATION, SAME PULSE DURATIONS =====")
  #Filters: double pulses, same wavelengths, same polarization
  PlottingDoublePulses = PlottingTable #get the full database
  PlottingDoublePulses = libDatabase.FilterDatabaseGreaterThan(PlottingDoublePulses, 1E0, libPlotting.IndexTable['Efield1']) #removes fields where Epeak < 1 V/m. #works
  PlottingDoublePulses = libDatabase.FilterDatabaseGreaterThan(PlottingDoublePulses, 1E0, libPlotting.IndexTable['Efield2']) #removes fields where Epeak2 < 1 V/m to select double pulses.
  PlottingDoublePulses = libDatabase.FilterDatabase(PlottingDoublePulses, float(0.), libPlotting.IndexTable['CrossPolarization']) #removes fields where Polarization1 != Polarization2
  PlottingDoublePulses = libDatabase.FilterDatabase(PlottingDoublePulses, PulseDuration, libPlotting.IndexTable['pulses2']) # Removing fields were pulseDuration1 != pulseDuration2

  
  #### 2D PLOTTING Density(Color1, Color2) for one given set of laser parameters. We produce one plot like this for equal fields. 
  
  ## Selecting fields which have a strength close to a number indicated in input. 
  PlottingDoublePulses = libDatabase.FilterDatabaseGreaterThan(PlottingDoublePulses, FieldStrength*(1.-Tolerance), libPlotting.IndexTable['Efield1'])
  PlottingDoublePulses = libDatabase.FilterDatabaseLowerThan(PlottingDoublePulses, FieldStrength*(1.+Tolerance), libPlotting.IndexTable['Efield1'])
  PlottingDoublePulses = libDatabase.FilterDatabaseGreaterThan(PlottingDoublePulses, FieldStrength*(1.-Tolerance), libPlotting.IndexTable['Efield2'])
  PlottingDoublePulses = libDatabase.FilterDatabaseLowerThan(PlottingDoublePulses, FieldStrength*(1.+Tolerance), libPlotting.IndexTable['Efield2'])

  
  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq, wavelength2_uniq, CEP2_uniq, Efield2_uniq, pulses2_uniq= libPlotting.CountTheNumberOfLaserParameters(PlottingDoublePulses)
  
  print(Header+"Total number of datasets: "+str(len(set(wavelength_uniq))*len(set(wavelength2_uniq))))
  
  #preparing the 2D plot
  DensityMax2D = np.empty((0,3))
  CurrentMax2D = np.empty((0,3))
  EnergyMax2D  = np.empty((0,3))
  Wpi2D        = np.empty((0,3))
  for Color1 in wavelength_uniq: #selecting color1
    PlottingDoublePulses_SelectedColor1 = libDatabase.FilterDatabase(PlottingDoublePulses, Color1, libPlotting.IndexTable['wavelengths1'])
    for Color2 in wavelength2_uniq: #selecting color2
      # Selecting Color 2
      PlottingDoublePulses_SelectedColor2 = libDatabase.FilterDatabase(PlottingDoublePulses_SelectedColor1, Color2, libPlotting.IndexTable['wavelengths2'])
      
      # Constructing the density table



      runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingDoublePulses_SelectedColor2)

      
      print("Color1: "+str(Color1)+" m, Color2: "+str(Color2)+" m.")
      if(len(Ne_final) >= 1):
        ToBeAdded        = np.array([wavelengths1[0], wavelengths2[0], Ne_final[0]])
        ToBeAddedCurrent = np.array([wavelengths1[0], wavelengths2[0], TotalCurrent_final[0]])
        ToBeAddedEnergy  = np.array([wavelengths1[0], wavelengths2[0], ElectronEnergy_max[0]])
        ToBeAddedWpi     = np.array([wavelengths1[0], wavelengths2[0], W_pi[0]])
        logger.debug(ToBeAdded, ToBeAddedCurrent, ToBeAddedEnergy, ToBeAddedWpi)
        DensityMax2D=np.vstack((DensityMax2D, ToBeAdded))
        CurrentMax2D=np.vstack((CurrentMax2D, ToBeAddedCurrent))
        EnergyMax2D =np.vstack((EnergyMax2D, ToBeAddedEnergy))
        Wpi2D       =np.vstack((Wpi2D, ToBeAddedWpi))
        del ToBeAdded, ToBeAddedCurrent, ToBeAddedEnergy, ToBeAddedWpi
        
  # Stack is built. Export, show and plot on a heat map. 
  logger.info("** Info: exporting the sub-plot2D: Ne_final(wavelength, Field), TotalCurrent_Final(), ElectronEnergy_max()...")
  logger.info("** Info: useful to visualize available data with GNUplot. ")
  filename         = "Bicolor-Field"+str(round(FieldStrength*1E-9))+"VperNm-NeFinal"
  filename_current = "Bicolor-Field"+str(round(FieldStrength*1E-9))+"VperNm-JeFinal"
  filename_energy  = "Bicolor-Field"+str(round(FieldStrength*1E-9))+"VperNm-Energy"
  filename_wpi     = "Bicolor-Field"+str(round(FieldStrength*1E-9))+"VperNm-Wpi"
  
  # Sorting the data before storing and plotting
  DensityMax2Dsorted = DensityMax2D[np.lexsort((DensityMax2D[:,2], DensityMax2D[:,1], DensityMax2D[:,0]))]
  CurrentMax2Dsorted = CurrentMax2D[np.lexsort((CurrentMax2D[:,2], CurrentMax2D[:,1], CurrentMax2D[:,0]))]
  EnergyMax2Dsorted  = EnergyMax2D[np.lexsort((EnergyMax2D[:,2], EnergyMax2D[:,1], EnergyMax2D[:,0]))]
  Wpi2Dsorted        = Wpi2D[np.lexsort((Wpi2D[:,2], Wpi2D[:,1], Wpi2D[:,0]))]
  #DensityMax2D.sort(axis=1) ##WARNING: THIS IS OUT OF CONTROL. DO NOT USE. 
  
  libDatabase.ExportToTxt(DensityMax2Dsorted, filename+".csv") #this is voluntarily unsorted to avoid damage
  libDatabase.ExportToTxt(CurrentMax2Dsorted, filename_current+".csv") #this is voluntarily unsorted to avoid damage
  libDatabase.ExportToTxt(EnergyMax2Dsorted,  filename_energy+".csv") #this is voluntarily unsorted to avoid damage
  libDatabase.ExportToTxt(Wpi2Dsorted,        filename_wpi+".csv") #this is voluntarily unsorted to avoid damage
  
  ## Lazy to build up a new function for that for now. 
  if(WhatToPlot == 1): 
      DensityMax2Dsorted = CurrentMax2Dsorted; 
      filename = filename_current
      QuantityTitle="Shift current $J_{end}(\lambda_1,\lambda_2)$"
  elif(WhatToPlot == 2): 
      DensityMax2Dsorted = EnergyMax2Dsorted
      filename = filename_energy
      QuantityTitle=r"Total energy $\xi_{end}(\lambda_1,\lambda_2)$"
  elif(WhatToPlot == 3): 
      DensityMax2Dsorted = Wpi2Dsorted
      filename = filename_wpi
      QuantityTitle=r"Excitation rate $w_{PI}(\lambda_1,\lambda_2)$"
  else: 
      QuantityTitle=r"Excited density $N_{exc}(\lambda_1,\lambda_2)$"
  
  #print "Sorted matrix?"
  #print DensityMax2Dsorted
  BrutX = 1E9*DensityMax2Dsorted[:,0]
  BrutY = 1E9*DensityMax2Dsorted[:,1]
  BrutZ = DensityMax2Dsorted[:,2]
  logger.debug(len(BrutX))
  logger.debug(len(BrutY))
  logger.debug(len(BrutZ))
  
  X=list(set(BrutX))
  Y=list(BrutY) ##BUG: we have degeneracy, investigate why and remove the problem. 
  Z=DensityMax2Dsorted
  #print X
  #print Y
  
  #print len(X)
  #print len(Y)
  YY, XX = np.meshgrid(X, Y)
  #print np.shape(XX)
  #print np.shape(YY)
  #print np.shape(DensityMax2Dsorted)
 
  #print DensityMax2Dsorted
  #print np.shape(XX)
  #DensityMax2DsortedShape = np.reshape(DensityMax2Dsorted, np.shape(XX))
  
  # Interpolation to smooth the picture and colorize it more softly
  # BUG: function interp2d does not work! 
  #MeshPrecision=50
  #DensityMax2D_interp = interp2d(BrutX, BrutY, BrutZ)
  #Xi = np.arange(np.min(BrutX), np.max(BrutX), (np.max(BrutX)-np.min(BrutX))/MeshPrecision)
  #Yi = np.arange(np.min(BrutY), np.max(BrutY), (np.max(BrutY)-np.min(BrutY))/MeshPrecision)
  #YYi, XXi = np.meshgrid(Xi, Yi)
  ##print XXi
  #BicolorNeFinal=DensityMax2D_interp(Xi, Yi)
  
  # If the size is big enough, let's make a plot.
  if(len(DensityMax2D) >= 2):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.title(QuantityTitle+', '+r'$E_{1,2}=$'+str(FieldStrength*1E-9)+' V/nm')
    plt.xlabel(r'Wavelength $\lambda_1$ (nm)')
    plt.ylabel(r'Wavelength $\lambda_2$ (nm)')
    #CS = plt.imshow(XX, YY, DensityMax2Dsorted) #, cmap=plt.cm.Blues)
    #CS = plt.contourf(XX, YY, DensityMax2Dsorted, cmap=plt.cm.Blues)
    #CS = ax.scatter(BrutX, BrutY, BrutZ) #, cmap=plt.cm.Blues) #WORKS
    #CS = ax.plot_wireframe(XX, YY, DensityMax2Dsorted) #, cmap=plt.cm.Blues) #FAILS: WHY do we have ground values ?! 
    #CS = ax.plot_surface(XX, YY, DensityMax2DsortedShape) #, cmap=plt.cm.Blues) #FAILS: WHY do we have ground values ?! 
    CS = ax.plot_trisurf(BrutX, BrutY, BrutZ, cmap=plt.cm.Blues) #WORKS
    
    #print np.shape(BrutX), np.shape(BrutY), np.shape(BrutZ)
    #print np.shape(Xi), np.shape(Yi)
    #print np.shape(XXi), np.shape(YYi)
    #print np.shape(BicolorNeFinal)
    
    #CS = ax.scatter(XXi, YYi, BicolorNeFinal) #, cmap=plt.cm.Blues) #WORKS
    
    #CS = ax.plot_trisurf(Xi, Yi, BicolorNeFinal, cmap=plt.cm.Spectral)
    
    #CS = ax.contour(XX, YY, Z, cmap=plt.cm.Blues) #WORKS
    #plt.colorbar(CS)
  
    plt.savefig(filename+'.eps')
    plt.savefig(filename+'.png')
    if(ShowFinalPlot==1):
      plt.show()
  
  del DensityMax2D, DensityMax2Dsorted
#}}}

def plotAdjustedCollisionFrequencies_wavelength(): #{{{
# Database was constructed by hand

    Efield = np.array(          [1E9,      2E9,       3E9,       4E9,       5E9,      1E10])
    CollTime_1200nm = np.array( [126E-15,  22.75E-15, 16.05e-15, 14.1e-15,  1.27e-14, 1.24e-14])
    CollTime_800nm = np.array(  [4.85e-15, 5.6E-15,   5.9E-15,   6.55E-15,  7.25e-15, 6.9e-15])
    CollTime_600nm = np.array(  [np.nan,   7.3E-15,   6.95E-15,  7.2E-15,   6.85e-15, 5.60e-15])
    CollTime_400nm = np.array(  [1.45e-15, 1.55e-15,  1.75E-15,  2.1E-15,   2.65e-15, 3.55e-15])

    # Select the runs where Collision Frequency is provided. 
    # Plot as function of the peak field amplitude. 

    plt.figure()
    plt.xlabel(r'Field peak amplitude $E_{\mathrm{peak}}$ (V/nm)')
    plt.ylabel(r'Estimated collision time $\nu^{-1}$ (fs)')
    plt.loglog(Efield*1E-9, 1E15*CollTime_400nm,  'b-+',  label=r"$\lambda=400$ nm")
    plt.loglog(Efield*1E-9, 1E15*CollTime_600nm,  'g-+',  label=r"$\lambda=600$ nm")
    plt.loglog(Efield*1E-9, 1E15*CollTime_800nm,  'k-+',  label=r"$\lambda=800$ nm")
    plt.loglog(Efield*1E-9, 1E15*CollTime_1200nm, 'r-+', label=r"$\lambda=1200$ nm") 
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig("CollisionFrequencies.eps")
    #plt.show()
    
#}}}                                                           

