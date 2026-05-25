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

## @module plotSinglePulse
# Plots for single-colour, single-pulse TDDFT runs.
# Functions: plotEffectOfRiseTime, plotEffectOfIntensity_wavelengths,
#            plotEffectOfIntensity_PulseDuration,
#            plotEffectOfIntensity_PulseDuration_Metals,
#            plotEffectOfIntensity_CEP, plotEffectOfCEP,
#            plotEffectOfFrozenXC

from octopus_slabs.Libs._plotCommon import *  # noqa: F401,F403

if EnablePeakFinder:
    from scipy.signal import find_peaks

def plotEffectOfRiseTime(PlottingTable, pulseWavelength=800e-9, pulseDuration=20e-15, Efield1=5E9, PulseShape=3, ShowFinalPlot=0, Egap_SI=2.56*e): #{{{
  logger.info("#====== Plot: SINGLE PULSE, EFFECT OF THE RISE TIME, 1 WAVELENGTH, 1 POLARIZATIONS: final density as function of RISETIME ======")
  PlottingSinglePulses = PlottingTable #get the full database
  PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, libPlotting.IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
  PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, libPlotting.IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0., libPlotting.IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, pulseWavelength, libPlotting.IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts. #TODO: this field should be making a loop in available wavelengths given in wavelengths1_uniq.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0, libPlotting.IndexTable['FreezeHXC']) #keeps only non-frozen XC.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, PulseShape, libPlotting.IndexTable['PulseType']) #keeps only non-frozen XC.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, pulseDuration, libPlotting.IndexTable['pulses1']) #removes fields where tau != pulseDuration

  PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, Efield1*1.01, libPlotting.IndexTable['Efield1'])
  PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, Efield1*0.99, libPlotting.IndexTable['Efield1'])

  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = libPlotting.CountTheNumberOfLaserParameters(PlottingSinglePulses)
  
  plt.figure(figsize=(SizeX,SizeY))
  ax1 = plt.subplot(111) #111 or plotA if using plotB
  plt.title(r'$n_{\mathrm{exc}}(RiseCycles)$, $\tau=$'+str(int(pulseDuration*1E15))+' fs.')
  # Number of sweeped parameters: 3 (Efield, Delay, CrossPolarization)
  #plt.scatter(libKeldysh.FieldToIntensity(Efield1), Ne_final, c=Delay, s=(30.+CrossPolarization*50.), label=r'$N_{exc}$[Field]')
  plot1 = [] #np.empty((0,0))
  for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
    PlottingSinglePulses_CEP = libDatabase.FilterDatabase(PlottingSinglePulses, CEP1_value, libPlotting.IndexTable['CEP1']) #select the local CEP
    #TODO: filter by field strength ! Maybe take a piece of code from tolerance-based selector in another function. #TODO: Use subfunctions otherwise we will get crazy. 
    ## prints the current state of selection and stores the values of sweeped parameters

    runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingTable)

    #plot_field, = ax1.plot(libKeldysh.FieldToIntensity(Efield1), Ne_final, '-', label=r'$N_{exc}^{final}$, $\phi_{CE}$='+str(round(CEP1_value*180./pi)))
    #lists = sorted(itertools.izip(*[libKeldysh.FieldToIntensity(Efield1), Ne_final]))
    lists = sorted(zip(*[RiseTime, Ne_max, Ne_final])) #NOTE: Ne_max makes more sense than Ne_final
    RiseTime_sorted, Ne_max_sorted, Ne_final_sorted = list(zip(*lists))
    #print Header+"Ready to plot: "
    #print "** Info: Sorting data..."
    #print Efield1_sorted, Ne_final_sorted
    plot_field, = ax1.plot(RiseTime_sorted, Ne_final_sorted, '-', label=r'$N_{exc}^{final}$, $\phi_{CE}$='+str(round(CEP1_value*180./pi)))
    logger.info(Header+"plotEffectOfIntensity_CEP(): adding CEP "+str(round(CEP1_value*180./pi))+" to the the plot.")
    logger.info(plot_field)
    #print Header+"** Debug info: "+str(np.shape(plot_field))
    plot1 += [plot_field] # does cumulate each CEP now. Type: Line2D, not np.array. 
    
    #ax1.legend(loc=4, borderaxespad=0.) #bbox_to_anchor=(1.05, 1), 
  labels1 = [l.get_label() for l in plot1]
  ax1.legend(plot1, labels1, loc='best', fancybox=True, framealpha=0)
  plt.xlabel(r'Rise cycles')
  plt.ylabel(r'Final density (m$^{-3}$)')
  #  plt.ylim((1E16,4.*5e28))
  #ax1.set_xscale('log')
  #ax1.set_yscale('log')
  #plot1_list = plot1.tolist()
    #print Efield1, CEP1_value
  
  ax12 = ax1.twinx()
  for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
    PlottingSinglePulses_CEP = libDatabase.FilterDatabase(PlottingSinglePulses, CEP1_value, libPlotting.IndexTable['CEP1']) #select the local CEP
    
    ## prints the current state of selection and stores the values of sweeped parameters
    runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingTable)
    
    #plot_field, = ax1.plot(libKeldysh.FieldToIntensity(Efield1), Ne_final, '-', label=r'$N_{exc}^{final}$, $\phi_{CE}$='+str(round(CEP1_value*180./pi)))
    #lists = sorted(itertools.izip(*[libKeldysh.FieldToIntensity(Efield1), Ne_final]))
    lists = sorted(zip(*[RiseTime, Ne_max, Ne_final, pulses1])) #NOTE: Ne_max makes more sense than Ne_final
    RiseTime_sorted, Ne_max_sorted, Ne_final_sorted, pulses1_sorted = list(zip(*lists))
    plt.plot(RiseTime_sorted, np.divide(Ne_final_sorted, pulses1_sorted), 'o') #, label=r'$N_{exc}/\tau$') #, I='+str(libKeldysh.FieldToIntensity(Efield1_value))+r' W/m$^2$')
  
  #ax12.set_yscale('log')
  ax12.set_ylabel(r'$w_{\mathrm{PI}}$ (m$^{-3}$.s$^{-1}$)')

  #ax2 = plt.subplot(plotB)
  #plt.title('Virtual excitation factor')
  
  #plt.plot(libKeldysh.FieldToIntensity(Efield1), irreversibility, 'ow', label='Virtual exc.')
  
  #plt.xlabel(r'Peak intensity (W/m$^2$)')
  #plt.ylabel('Virtual exc. factor')
  #plt.ylim((0,1))
  #ax2.set_xscale('log')
  ##adding the Keldysh adiabatic factor
  #ax22 = ax2.twinx()
  
  #plt.plot(libKeldysh.FieldToIntensity(Efield1), libKeldysh.gammaKeldysh(2.56*e, 0.228, Efield1, wavelengths1), 'or', label='Adiab. fact.') #BUG: 2.56 must not be specified by hand here! 

  #plt.ylabel('Adiabadicity factor')
  ##plt.ylim((0,1))
  #ax22.set_yscale('log')

  #ax2.set_yscale("log")
  #ax2.legend(loc=2, borderaxespad=0.)
  #ax1.legend(bbox_to_anchor=(1.15, 1), loc=2, borderaxespad=0.)
  #ax22.legend(loc=1, borderaxespad=0.)

  # Saving the plot
  ExportName = "EffectOfRiseTime"
  plt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.eps')
  plt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.png')
  if(ShowFinalPlot==1):
    plt.show()

  del ax1, PlottingSinglePulses, wavelength_uniq, wavelengths1
#}}}

## Plots Nexc(intensity; single tau, 1 curve per wavelength) | Wpi(intensity)
def plotEffectOfIntensity_wavelengths(PlottingTable, pulseDuration=20e-15, CEP=0., RestrictWavelength1=False, wavelength_value=800e-9, PulseType=3, RiseTime=2, ShowFinalPlot=0, Egap_SI = 2.56*e): #{{{
  PlotVirtualExcitation = False
  print("#====== Plot: SINGLE PULSE, VARIOUS WAVELENGTHS, SAME PULSE DURATION, SAME POLARIZATIONS: final density as function of field ====== ")
  #Filter: double pulses, same wavelengths, same polarization
  PlottingSinglePulses = PlottingTable #get the full database
  PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, libPlotting.IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
  PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, libPlotting.IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0., libPlotting.IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, CEP, libPlotting.IndexTable['CEP1']) #removes fields where CEP != 0. #NOTE: Can be done for Si-LDA.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0., libPlotting.IndexTable['FreezeHXC']) #removes fields where FreezeHXC != 0.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0., libPlotting.IndexTable['FieldInduced']) #removes fields where FreezeHXC != 0.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 1., libPlotting.IndexTable['TheoryLevel']) #removes runs where TheoryLevel != dft.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, PulseType, libPlotting.IndexTable['PulseType']) #keeps only pulses of type StH (3)
  #PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, RiseTime, libPlotting.IndexTable['RiseTime']) #keeps only pulses of type STH RiseTime
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, pulseDuration, libPlotting.IndexTable['pulses1']) #keeps only pulses of selected duration in argument
  if(RiseTime == 2 or RiseTime==0): #select by RiseTime only if it is defined. I include only 2 and 0 here. 
      PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, RiseTime, libPlotting.IndexTable['RiseTime'])
  
  ## prints the current state of selection
  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq, wavelength2_uniq,CEP2_uniq, \
  Efield2_uniq,pulses2_uniq = libPlotting.CountTheNumberOfLaserParameters(PlottingSinglePulses)
  if(PlotVirtualExcitation==True): 
      plotA=211; plotB=212
  else: 
      plotA=111; plotB=111
  # Simple pulse, single color
  plt.figure()
  #if(len(wavelength_uniq) == 1): #single color simulations
  # Complexity level: 1
  plot1=[]
  ax1 = plt.subplot(plotA)
  
  plt.title('Density as function of peak intensity')
  # Number of sweeped parameters: 3 (Efield, Delay, CrossPolarization)
  #plt.scatter(libKeldysh.FieldToIntensity(Efield1), Ne_final, c=Delay, s=(30.+CrossPolarization*50.), label=r'$N_{exc}$[Field]')
  for wavelength_value in wavelength_uniq: #{{{ #TODO: wavelength_uniq should be sorted (visible in the caption)
    #if(RestrictWavelength1==False): 
      # select the corresponding wavelength in results
    PlottingSinglePulses_wavelength = libDatabase.FilterDatabase(PlottingSinglePulses, wavelength_value, libPlotting.IndexTable['wavelengths1']) #keeps all wavelengths wavelenghts.
    #else: #selection only 1 wavelength #BUG: DOES NOT WORK, and even attributes to other wavelength what was done for 800 nm
      #PlottingSinglePulses_wavelength = libDatabase.FilterDatabase(PlottingSinglePulses, 800e-9, IndexTable['wavelengths1']) #keeps only 800 nm

    runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingTable)
    lists = sorted(zip(*[wavelengths1,libKeldysh.FieldToIntensity(Efield1), Ne_final, pulses1]))

    wavelengths_sorted, Intensity_sorted, Ne_final_sorted, pulses1_sorted = list(zip(*lists))
    plotNe, = plt.plot(Intensity_sorted, Ne_final_sorted, 'o', label=r'$N_{exc}$, $\lambda=$'+str(round(wavelength_value*1E9,0))+' nm')
  #}}}
  #ax1.legend(loc='best', fancybox=True, framealpha=0) #bbox_to_anchor=(1.05, 1), 
  if(not PlotVirtualExcitation): 
    plt.xlabel(r'Peak intensity (W/cm$^2$)')
  plt.ylabel(r'Final exc. density (m$^{-3}$)')
#  plt.ylim((1E16,4.*5e28))
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  
  ax12 = ax1.twinx()
  for wavelength_value in wavelength_uniq: #{{{
    # select the corresponding wavelength in results
    if(RestrictWavelength1==False): 
      PlottingSinglePulses_wavelength = libDatabase.FilterDatabase(PlottingSinglePulses, 800e-9, libPlotting.IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts.
    else:
      PlottingSinglePulses_wavelength = libDatabase.FilterDatabase(PlottingSinglePulses, wavelength_value, libPlotting.IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts.

    runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingTable)

    lists = sorted(zip(*[wavelengths1,libKeldysh.FieldToIntensity(Efield1), Ne_final, pulses1]))
    wavelengths_sorted, Intensity_sorted, Ne_final_sorted, pulses1_sorted = list(zip(*lists))

    plt.plot(Intensity_sorted, np.divide(Ne_final_sorted, (np.multiply(2.,pulses1_sorted))), 'o') #, label=r'$N_{exc}/\tau$, $\lambda=$'+str(round(wavelength_value*1E9))+' nm')
  #}}}
  ax12.set_yscale('log')
  plt.ylabel('$w_{\mathrm{PI}}$ (m$^{-3}$.s$^{-1}$)')
  
  # Second panel of the figure
  if(PlotVirtualExcitation):
    ax2 = plt.subplot(plotB)
    #plt.title('rVirtual exc. $\zeta$')
    plt.plot(libKeldysh.FieldToIntensity(Efield1), irreversibility, 'ob', label='Virtual exc.')
    plt.xlabel(r'Peak intensity (W/cm$^2$)')
    plt.ylabel(r'Virtual exc. $\zeta$')
    plt.ylim((0,1))
    ax2.set_xscale('log')
    ax22 = ax2.twinx()
    for wavelength_value in wavelength_uniq:
      # select the corresponding wavelength in results
      PlottingSinglePulses_wavelength = libDatabase.FilterDatabase(PlottingSinglePulses, wavelength_value, libPlotting.IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts.

      runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingSinglePulses_wavelength)
      plt.plot(libKeldysh.FieldToIntensity(Efield1), libKeldysh.gammaKeldysh(Egap_SI, meff, Efield1, wavelengths1), 'or', label=r'Adiabaticity, $\lambda=$'+str(round(wavelength_value*1E9, 0))+' nm')

    plt.ylabel(r'Adiabadicity $\gamma$')
    #plt.ylim((0,1))
    ax22.set_yscale('log')
      
    # ax2.set_yscale("log")
    ax2.legend(loc='best', fancybox=True, framealpha=0)
    # ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax22.legend(loc='best', fancybox=True, framealpha=0)
  ax12.legend(loc='best', fancybox=True, framealpha=0)
  
  # Saving the plot
  ExportName = "EffectOfIntensity-"+str(int(pulseDuration*1E15))+"fs"
  plt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.eps')
  plt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.png')
  logger.info(Header+"** Saved: Metaplot-SinglePulse-SingleColor-"+ExportName+".eps|png.")
  if(ShowFinalPlot==1):
    plt.show()

  del ax1, ax12, PlottingSinglePulses, wavelength_uniq, CEP1_uniq
  if(PlotVirtualExcitation): 
    del ax2, ax22
#}}}

## Plots (N_exc | Wpi | Energy) (peak intensity), (1 curve per pulse duration)
# Nexc(I) is used to extract 1 to 7 photon absorption rates, given in caption. 
# Wpi is compared to Keldysh analytical solid theory, and numerically integrated Keldysh if available. # Wpi is also compared to tunneling limit case. 
# ElectronEnergy (intensity) is given. 
def plotEffectOfIntensity_PulseDuration(PlottingTable, pulseWavelength=800e-9, CEP0=0., PulseType=3, RiseCycles=2, PlotWpiAxis=False, ShowFinalPlot=0, Egap_SI=2.56*e, Egap_SI_indirect=0.51*e, OpticalIndex=1.0, MoveIon=0, Functional="LDA"): #{{{
    UseNexcAfterPulse = False #instead of density at end of simulation, one capture the ones very close to the end of the pulse.
    PlotEnergy        = False #PRL: False. #PRB: True.
    # EnableKeldysh     = True
    # EnableStarkEffect =False #Stark effect on tunneling? Disabled in original Keldysh, Gruzdev and Kaiser theories. Now that we know Stark is included in Keldysh, it does not make sense to enable it for tunneling limit case (Stark ionization potential is time-average of the dynamically triggered electron motion).
    PlotSaturation    =False
    ShowPulseDurations=True #True #: show the pulse duration in caption. False: hides it (for saving visual space)
    WithOpticalIndex  =False  #True for solids. False for atomic Keldysh. #NOTE: this curious feature may originate from the usage of atomic units in the effective mass approximation (just an hypothesis).
    SidePanels        =True #True: [a] | [b]. False: a above b. 
    ShowSigmaRates    = 2 #Cosmetics: 0: remove all fittings, 1: show n=..., 2: show sigma values
    EnableInsetMagnification=False #this is defined as function of wavelength below
    PlotVirtualExcitation = False #Not advised! This is not a real observable.

    global meff_Waroquiers_LDA
    Ntotal=5.0E28
    meff_Waroquiers_LDA = 0.2226 #0.2226 #Waroquiers, D. et al. Physical Review B, 2013, 87, 075121.
    meff = meff_Waroquiers_LDA #0.14 #used to adjust the w_PI model
    # meff_Tunnel = meff_Waroquiers_LDA #fits well if using INDIRECT band gap.

    # This disables the work made on the optical index. Important to show in the Suppl. Inf.
    if(not WithOpticalIndex):
        OpticalIndex=1.
        FileName_opt="-NoOpticalIndex"
    else:
        FileName_opt=""

    x_min  = {407.58E-9: (1E9, 1E15), 483e-9: (1E9, 1E15), 800e-9: (1E9, 1E15), 1030e-9: (1E9, 1E15), 1600e-9: (1E10, 1E15), 3200e-9: (1E10, 1E15)}

    print("====== Plot: SINGLE PULSE, SAME WAVELENGTHS, SAME POLARIZATIONS: final density as function of field and pulse duration ======")
    #Filter: double pulses, same wavelengths, same polarization
    PlottingSinglePulses = PlottingTable #get the full database
    PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, libPlotting.IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
    PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, libPlotting.IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
    PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, libPlotting.IndexFunctional[Functional], libPlotting.IndexTable['Functional']) #removes runs where functional is not the chosen functional
    PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0., libPlotting.IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0.
    PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, CEP0, libPlotting.IndexTable['CEP1']) #removes fields where CEP != 0. #NOTE: Can be done for Si-LDA.
    #PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, CEP0, libPlotting.IndexTable['CEP1'])
    PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0., libPlotting.IndexTable['FreezeHXC']) #removes fields where FreezeHXC != 0.
    PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0., libPlotting.IndexTable['FieldInduced']) #removes fields where FreezeHXC != 0.
    PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 1., libPlotting.IndexTable['TheoryLevel']) #removes runs where TheoryLevel != dft.
    PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, PulseType, libPlotting.IndexTable['PulseType']) #keeps only pulses of type StH (3)
    PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, RiseCycles, libPlotting.IndexTable['RiseTime']) #keeps only pulses of a selected RiseCycles
    PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, pulseWavelength, libPlotting.IndexTable['wavelengths1']) #keeps only pulses of type STH (2)
    PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, MoveIon, libPlotting.IndexTable['MoveIon']) #filters the enabled or disabled movement of ions.
    wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq, wavelength2_uniq, CEP2_uniq, Efield2_uniq, pulses2_uniq = libPlotting.CountTheNumberOfLaserParameters(
        PlottingSinglePulses)

    #from collections import OrderedDict

    #linestyles_dict = OrderedDict(
    #[('solid',               (0, ())),
    #('loosely dotted',      (0, (1, 10))),
    #('dotted',              (0, (1, 5))),
    #('densely dotted',      (0, (1, 1))),

    #('loosely dashed',      (0, (5, 10))),
    #('dashed',              (0, (5, 5))),
    #('densely dashed',      (0, (5, 1))),

    #('loosely dashdotted',  (0, (3, 10, 1, 10))),
    #('dashdotted',          (0, (3, 5, 1, 5))),
    #('densely dashdotted',  (0, (3, 1, 1, 1))),

    #('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
    #('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
    #('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))])

    ## We simply select the pulse duration larger than 21 fs (to not overload pics).
    if(pulseWavelength==3200e-9): #NOTE: enable if you wish to remove 20 fs results for 3200 nm.
        #PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, 15e-15, libPlotting.IndexTable['pulses1'])
        try:
            PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, 9e-15, libPlotting.IndexTable['pulses1'])
        except:
            logger.error(Header + "** Warning: no data for " + str(pulseWavelength))
            os._exit(1)
        #EnableInsetMagnification #Show magnification only for matching results
    if(pulseWavelength==1600e-9): #NOTE: enable if you wish to remove 20 fs results for 1600 nm.
        #PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, 15e-15, libPlotting.IndexTable['pulses1'])
        try:
            PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, 9e-15, libPlotting.IndexTable['pulses1'])
        except:
            logger.error(Header + "** Warning: no data for " + str(pulseWavelength))
            os._exit(1)
        #EnableInsetMagnification=False
    if (pulseWavelength == 1030e-9):  # NOTE: enable if you wish to remove 20 fs results for 1600 nm.
        # PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, 15e-15, libPlotting.IndexTable['pulses1'])
        try:
            PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, 9e-15,
                                                                         libPlotting.IndexTable['pulses1'])
            PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, 21e-15,
                                                                         libPlotting.IndexTable['pulses1'])
        except:
            logger.error(Header + "** Warning: no data for " + str(pulseWavelength))
            os._exit(1)
    if(pulseWavelength==800e-9): #NOTE: enable if you wish to remove 10 fs results for 800 nm.
        #PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, 15e-15, libPlotting.IndexTable['pulses1'])
        try:
            PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, 9e-15, libPlotting.IndexTable['pulses1'])
        except:
            logger.error(Header + "** Warning: no data for " + str(pulseWavelength))
            os._exit(1)
        #EnableInsetMagnification=False
    if(pulseWavelength==483e-9): #NOTE: enable if you wish to remove 10 fs results for 800 nm.
        #PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, 15e-15, libPlotting.IndexTable['pulses1'])
        try:
            PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, 9e-15, libPlotting.IndexTable['pulses1'])
        except:
            logger.error(Header + "** Warning: no data for " + str(pulseWavelength))
            os._exit(1)
        #EnableInsetMagnification=False

    ## NOTE: For clarity, I will filter results that are supplementary for the paper's figures.
    ## Let's remove the results of 1600 nm with 20 fs pulse duration. We can keep 30 fs.

    ## prints the current state of selection
    wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = libPlotting.CountTheNumberOfLaserParameters(PlottingSinglePulses) #TODO: use PANDAS here

    # Simple pulse, single color
    #plotA = 311; plotB = 312; plotC = 313 #Horizontal
    if(not PlotEnergy):
        if(SidePanels):
            plotA=121; plotB=122;#2 columns figure
            ncols=1
        else:
            plotA=211; plotB=212;#2 rows figure
            ncols=2
            #plt.figure(figsize=(SizeY, SizeX/3.*2.))
        plt.figure(figsize=(SizeX/3.*2.,SizeY))
    else:
        plotA = 131; plotB = 132; plotC = 133 #Vertical
        plt.figure(figsize=(SizeX,SizeY))

    #plotA = 221; plotB = 224; plotC = 223; plotD = 222 # 2x2 for huge caption of ax1

    # EnergyMin = 0e0 #-1E-7
    NumberOfPhotons_round = 99 #we need to initialize it here

    ax1 = plt.subplot(plotA)

    plot1 = []
    Symbols={10e-15: '+',      1.03071305416229e-14: '+', 10.674e-15: '+', 20e-15: 'o', 30e-15: '^'}
    Colors={10e-15: 'g',       1.03071305416229e-14: 'g', 10.674e-15: 'g', 20e-15: 'r', 30e-15: 'b'}
    LineWidth={10e-15: '10.',  1.03071305416229e-14: '10.', 10.674e-15: '10.', 20e-15: '0.1', 30e-15: '2.'}

    linestyle_tuple = {
        'solid':                (0, ()),
        'dotted':               (0, (1,1)),
        'loosely dotted':        (0, (1, 10)),
        'dotted':                (0, (1, 1)),
        'densely dotted':        (0, (1, 1)),
        'loosely dashed':        (0, (5, 10)),
        'dashed':                (0, (5, 5)),
        'densely dashed':        (0, (5, 1)),
        'loosely dashdotted':    (0, (3, 10, 1, 10)),
        'dashdotted':            (0, (3, 5, 1, 5)),
        'densely dashdotted':    (0, (3, 1, 1, 1)),
        'dashdotdotted':         (0, (3, 5, 1, 5, 1, 5)),
        'loosely dashdotdotted': (0, (3, 10, 1, 10, 1, 10)),
        'densely dashdotdotted': (0, (3, 1, 1, 1, 1, 1))
    }
    ListOfInterestingGammas = [0.10, 1.0, 10.]

    for pulse_value in pulses_uniq: #{{{
        # select the corresponding wavelength in results
        PlottingSinglePulses_pulses = libDatabase.FilterDatabase(PlottingSinglePulses, pulse_value, libPlotting.IndexTable['pulses1'])
        PlottingSinglePulses_pulses = libDatabase.FilterDatabase(PlottingSinglePulses_pulses, pulse_value, libPlotting.IndexTable['pulses2']) #keeps only pulses of the same duration.

        runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIonRun, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingSinglePulses_pulses) #TODO: use PANDAS HERE

        logger.info(Header+"** Info: saving the data in CSV file for fitting. ")
        filename = "EffectOfPulseIntensity-Duration-"+str(pulseWavelength*1E9)+"nm-tau"+str(int(pulse_value*1E15))+"fs"+FileName_opt
        logger.info(Header+"** Info: saving in "+filename+".csv")
        libDatabase.ExportToTxt(PlottingSinglePulses_pulses, filename+".csv")
        if(UseNexcAfterPulse):
            del Ne_final
            Ne_final = Ne_afterpulse #NOTE: dirty and valid only until next call of ExtractData
        lists = sorted(zip(*[pulses1, Efield1, libKeldysh.FieldToIntensity(Efield1), Ne_final]))
        pulses1_sorted, Efield1_sorted, Intensity_sorted, Ne_final_sorted = list(zip(*lists))
        #print Header+"Ready to plot: "

        if(ShowPulseDurations):
            tau_text=r'$\tau=$'+str(int(pulse_value*1E15))+' fs'
        else:
            tau_text=""

        plot_field, = ax1.plot(np.multiply(xscale,Intensity_sorted), Ne_final_sorted,
                               Colors[pulse_value]+Symbols[pulse_value], linewidth=LineWidth[pulse_value], label=r'$n_{\mathrm{exc}}^{\mathrm{TDDFT}}$'+tau_text)
        plot1 += [plot_field]
        #plt.xlim(x_min[pulseWavelength])

        #==== CURVE FITTING SHOULD COME HERE === as it is for each pulse duration.
        logger.info(Header+"** Info: Computation of multi-photon excitation")

        ### Multiphotonic absorption law (used for fitting td-dft results)
        # d e f MultiPhotonAbsorption(intensity, sigma_n):

        #n_photons=2 # is defined by inheriting its value
        #omega = 2.*np.pi*c/wavelengths1[0] #how come this value is taken outside, but not my n_photon quantity?
        ##n_photons #gets the number of photons in upper levels
        #return sigma_n * intensity**n_photons / (n_photons*hbar*omega) * pulses1_sorted[0]

        ## 1-photon absorption fitting function
        def OnePhotonAbsorption(intensity, sigma_n):
            #{{{
            n=1; #wavelengths1=800e-9; #is captured by parental legacy
            omega = 2.*np.pi*c/wavelengths1[0]
            #nexc_mpi = MultiPhotonAbsorption(intensity, sigma_n, n, wavelengths1[0], pulses1_sorted[0])
            nexc_mpi = sigma_n * intensity**n / (n*hbar*omega) * pulses1_sorted[0]
            return nexc_mpi
        #}}}

        ## 2-photon absorption fitting function
        def TwoPhotonAbsorption(intensity, sigma_n):
            #{{{
            n=2; #wavelengths1=800e-9; #is captured by parental legacy
            omega = 2.*np.pi*c/wavelengths1[0]
            #nexc_mpi = MultiPhotonAbsorption(intensity, sigma_n, n, wavelengths1[0], pulses1_sorted[0])
            nexc_mpi = sigma_n * intensity**n / (n*hbar*omega) * pulses1_sorted[0]
            return nexc_mpi
        #}}}

        ## 3-photon absorption fitting function
        def ThreePhotonAbsorption(intensity, sigma_n):
            #{{{
            n=3; #wavelength=800e-9; tau_FWHM=20e-15;
            omega = 2.*np.pi*c/wavelengths1[0]
            tau_FWHM = pulses1_sorted[0]
            #nexc_mpi = MultiPhotonAbsorption(intensity, sigma_n, n, wavelength, tau_FWHM)
            nexc_mpi = sigma_n * intensity**n / (n*hbar*omega) * tau_FWHM
            return nexc_mpi
        #}}}

        ## 4-photon absorption fitting function
        def FourPhotonAbsorption(intensity, sigma_n):
            #{{{
            n=4; #wavelength=800e-9; tau_FWHM=20e-15;
            omega = 2.*np.pi*c/wavelengths1[0]
            tau_FWHM = pulses1_sorted[0]
            #nexc_mpi = MultiPhotonAbsorption(intensity, sigma_n, n, wavelength, tau_FWHM)
            nexc_mpi = sigma_n * intensity**n / (n*hbar*omega) * tau_FWHM
            return nexc_mpi
        #}}}

        ## 5-photon absorption fitting function
        def FivePhotonAbsorption(intensity, sigma_n):
            #{{{
            n=5; #wavelength=800e-9; tau_FWHM=20e-15;
            omega = 2.*np.pi*c/wavelengths1[0]
            tau_FWHM = pulses1_sorted[0]
            #nexc_mpi = MultiPhotonAbsorption(intensity, sigma_n, n, wavelength, tau_FWHM)
            nexc_mpi = sigma_n * intensity**n / (n*hbar*omega) * tau_FWHM
            return nexc_mpi
        #}}}

        ## 6-photon absorption fitting function
        def SixPhotonAbsorption(intensity, sigma_n):
            #{{{
            n=6; #wavelength=800e-9; tau_FWHM=20e-15;
            omega = 2.*np.pi*c/wavelengths1[0]
            tau_FWHM = pulses1_sorted[0]
            #nexc_mpi = MultiPhotonAbsorption(intensity, sigma_n, n, wavelength, tau_FWHM)
            nexc_mpi = sigma_n * intensity**n / (n*hbar*omega) * tau_FWHM
            return nexc_mpi
        #}}}

        ## 7-photon absorption fitting function
        def SevenPhotonAbsorption(intensity, sigma_n):
            #{{{
            n=7; #wavelength=800e-9; tau_FWHM=20e-15;
            omega = 2.*np.pi*c/wavelengths1[0]
            tau_FWHM = pulses1_sorted[0]
            #nexc_mpi = MultiPhotonAbsorption(intensity, sigma_n, n, wavelength, tau_FWHM)
            nexc_mpi = sigma_n * intensity**n / (n*hbar*omega) * tau_FWHM
            return nexc_mpi
        #}}}

        ### 2-photon and 3-photon absorption fitting function
        #def TwoAndThreePhotonAbsorption(intensity, sigma_2, sigma_3):
        ##{{{
        #n=2; #wavelengths1=800e-9; #is captured by parental legacy
        #omega = 2.*np.pi*c/wavelengths1[0]
        ##nexc_mpi = MultiPhotonAbsorption(intensity, sigma_n, n, wavelengths1[0], pulses1_sorted[0])
        #nexc_mpi = sigma_2 * intensity**2 / (2*hbar*omega) * pulses1_sorted[0] + sigma_3 * intensity**3 / (3*hbar*omega) * pulses1_sorted[0]
        #return nexc_mpi
        ##}}}

        #LineType_MPIfitArray={1: '--', 2: ':', 3: '-.', 4: ':', 5: '-'} #, 6: ':', 7: '--'}
        #LineType_MPIfitArray={1: linestyles_dict['solid'], 2: linestyles_dict['loosely dotted'], 3: linestyles_dict['dotted'], 4: linestyles_dict['densely dotted'], 5: linestyles_dict['dashed']}
        LineType_MPIfitArray={1: 'solid', 2: 'dotted', 3: 'dashed', 4: 'dashdotted', 5: 'loosely dashed', 6: 'densely dashed', 7: 'densely dashdotted'} #see https://matplotlib.org/3.1.0/gallery/lines_bars_and_markers/linestyles.html

        ## Filters a part of the data and fit using the given function
        def FitOnePhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, n_photons, OnePhotonIntensityRegime_min, OnePhotonIntensityRegime_max, plotMe=True): #TODO: pass function name in argument
            #{{{
            logger.info(Header+"** Info: Filtering the fitting region...")
            Intensity_sorted_np = np.array(Intensity_sorted) #filtering works for NumPy array, not for lists.
            #condition = (Intensity_sorted_np > OnePhotonIntensityRegime_min) and (Intensity_sorted_np < OnePhotonIntensityRegime_max)
            Intensity_sorted_mask = Intensity_sorted_np[(Intensity_sorted_np > OnePhotonIntensityRegime_min) & (Intensity_sorted_np < OnePhotonIntensityRegime_max)]

            Ne_final_sorted_np  = np.array(Ne_final_sorted)
            Ne_final_sorted_mask  = Ne_final_sorted_np[(Intensity_sorted_np > OnePhotonIntensityRegime_min) & (Intensity_sorted_np < OnePhotonIntensityRegime_max)]
            logger.info(Header+""+str(Intensity_sorted_mask)+""+str(Ne_final_sorted_mask)) #Works well

            OnePhotonAbsorption_params = curve_fit(OnePhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
            #TwoPhotonAbsorption_params = curve_fit(MultiPhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
            logger.info(Header+"** Ionization parameters for tau="+str(int(1e15*pulse_value))+" fs: "+str(OnePhotonAbsorption_params[0]))

            logger.info(Header+"** Info: adding the fitted curve to the plot")
            Nexc_1PI = OnePhotonAbsorption(Intensity_sorted_mask, OnePhotonAbsorption_params[0])
            #print Nexc_2PI
            if(ShowPulseDurations):
                tau_text=r', $\tau=$'+str(int(pulse_value*1E15))+' fs'
            else:
                tau_text=""
            if(ShowSigmaRates==2):
                sigma_txt=r"$\sigma_1=$"+str('%.4g'%OnePhotonAbsorption_params[0][0])+tau_text
            elif(ShowSigmaRates==1):
                sigma_txt=r"$n=1$"
            else:
                sigma_txt=""
            if(plotMe):
                plot_field, = ax1.plot(np.multiply(xscale,Intensity_sorted_mask), Nexc_1PI, Colors[pulse_value], linestyle=linestyle_tuple[LineType_MPIfitArray[1]], label=sigma_txt) #NOTE add value of sigma_1 here
                #ax1.fill_betweenx(Nexc_1PI, OnePhotonIntensityRegime_min, OnePhotonIntensityRegime_max, color='blue', alpha=0.5)
            else:
                plot_field = []
            del Nexc_1PI, OnePhotonAbsorption_params, Intensity_sorted_mask
            return plot_field
        #}}}

        ## Filters a part of the data and fit using the given function
        def FitTwoPhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, n_photons, TwoPhotonIntensityRegime_min, TwoPhotonIntensityRegime_max, plotMe=True): #TODO: pass function name in argument --> Pointer?
            #{{{
            logger.info(Header+"** Info: Filtering the fitting region...")
            Intensity_sorted_np = np.array(Intensity_sorted) #filtering works for NumPy array, not for lists.
            #condition = (Intensity_sorted_np > TwoPhotonIntensityRegime_min) and (Intensity_sorted_np < TwoPhotonIntensityRegime_max)
            Intensity_sorted_mask = Intensity_sorted_np[(Intensity_sorted_np > TwoPhotonIntensityRegime_min) & (Intensity_sorted_np < TwoPhotonIntensityRegime_max)]

            Ne_final_sorted_np  = np.array(Ne_final_sorted)
            Ne_final_sorted_mask  = Ne_final_sorted_np[(Intensity_sorted_np > TwoPhotonIntensityRegime_min) & (Intensity_sorted_np < TwoPhotonIntensityRegime_max)]
            # logger.info(Intensity_sorted_mask, Ne_final_sorted_mask) #Works well
            logger.info(Header + "" + str(Intensity_sorted_mask) + "" + str(Ne_final_sorted_mask))  # Works well
            TwoPhotonAbsorption_params = curve_fit(TwoPhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
            #TwoPhotonAbsorption_params = curve_fit(MultiPhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
            logger.info(Header+"** Ionization parameters for tau="+str(int(1e15*pulse_value))+" fs: "+str(TwoPhotonAbsorption_params[0]))

            logger.info(Header+"** Info: adding the fitted curve to the plot")
            Nexc_2PI = TwoPhotonAbsorption(Intensity_sorted_mask, TwoPhotonAbsorption_params[0])
            #print Nexc_2PI
            if(ShowPulseDurations):
                tau_text=r', $\tau=$'+str(int(pulse_value*1E15))+' fs'
            else:
                tau_text=""

            if(ShowSigmaRates==2):
                sigma_txt=r"$\sigma_2=$"+str('%.4g'%TwoPhotonAbsorption_params[0][0])+tau_text
            elif(ShowSigmaRates==1):
                sigma_txt=r"$n=2$"
            else:
                sigma_txt=""
            if(plotMe):
                plot_field, = ax1.plot(np.multiply(xscale,Intensity_sorted_mask), Nexc_2PI, Colors[pulse_value], linestyle=linestyle_tuple[LineType_MPIfitArray[2]], label=sigma_txt) #NOTE add value of sigma_2 here
            else:
                plot_field = []
            del Nexc_2PI, TwoPhotonAbsorption_params, Intensity_sorted_mask
            return plot_field
        #}}}


        ## ==== FITTING THE 3-photon ionization plot
        def FitThreePhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, n_photons, ThreePhotonIntensityRegime_min, ThreePhotonIntensityRegime_max, plotMe=True): #TODO: pass function name in argument instead of copying this function... :-(
            #{{{
            logger.info(Header+"** Info: Filtering the fitting region...")
            Intensity_sorted_np = np.array(Intensity_sorted) #filtering works for NumPy array, not for lists.
            #condition = (Intensity_sorted_np > TwoPhotonIntensityRegime_min) and (Intensity_sorted_np < TwoPhotonIntensityRegime_max)
            Intensity_sorted_mask = Intensity_sorted_np[(Intensity_sorted_np > ThreePhotonIntensityRegime_min) & (Intensity_sorted_np < ThreePhotonIntensityRegime_max)]

            Ne_final_sorted_np  = np.array(Ne_final_sorted)
            Ne_final_sorted_mask  = Ne_final_sorted_np[(Intensity_sorted_np > ThreePhotonIntensityRegime_min) & (Intensity_sorted_np < ThreePhotonIntensityRegime_max)]
            # logger.info(Intensity_sorted_mask, Ne_final_sorted_mask) #Works well
            logger.info(Header + "" + str(Intensity_sorted_mask) + "" + str(Ne_final_sorted_mask))  # Works well

            ThreePhotonAbsorption_params = curve_fit(ThreePhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
            #TwoPhotonAbsorption_params = curve_fit(MultiPhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
            logger.info(Header+"** Ionization parameters for tau="+str(int(1e15*pulse_value))+" fs: "+str(ThreePhotonAbsorption_params[0]))

            if(ShowPulseDurations):
                tau_text=r', $\tau=$'+str(int(pulse_value*1E15))+' fs'
            else:
                tau_text=""

            if(ShowSigmaRates==2):
                sigma_txt=r"$\sigma_3=$"+str('%.4g'%ThreePhotonAbsorption_params[0][0])+tau_text
            elif(ShowSigmaRates==1):
                sigma_txt=r"$n=3$"
            else:
                sigma_txt=""

            logger.info(Header+"** Info: adding the fitted curve to the plot")
            Nexc_3PI = ThreePhotonAbsorption(Intensity_sorted_mask, ThreePhotonAbsorption_params[0])
            #print Nexc_2PI
            if(plotMe):
                plot_field, = ax1.plot(np.multiply(xscale,Intensity_sorted_mask), Nexc_3PI, Colors[pulse_value], linestyle=linestyle_tuple[LineType_MPIfitArray[3]], label=sigma_txt) #NOTE adds value of sigma_3 here
            else:
                plot_field = []
            del Nexc_3PI, ThreePhotonAbsorption_params, Intensity_sorted_mask
            return plot_field
        #}}}

        def FitFourPhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, n_photons, FourPhotonIntensityRegime_min, FourPhotonIntensityRegime_max, plotMe=True): #TODO: pass function name in argument instead of copying this function... :-(
            #{{{
            logger.info(Header+"** Info: Filtering the fitting region...")
            Intensity_sorted_np = np.array(Intensity_sorted) #filtering works for NumPy array, not for lists.
            #condition = (Intensity_sorted_np > TwoPhotonIntensityRegime_min) and (Intensity_sorted_np < TwoPhotonIntensityRegime_max)
            Intensity_sorted_mask = Intensity_sorted_np[(Intensity_sorted_np > FourPhotonIntensityRegime_min) & (Intensity_sorted_np < FourPhotonIntensityRegime_max)]

            Ne_final_sorted_np  = np.array(Ne_final_sorted)
            Ne_final_sorted_mask  = Ne_final_sorted_np[(Intensity_sorted_np > FourPhotonIntensityRegime_min) & (Intensity_sorted_np < FourPhotonIntensityRegime_max)]
            logger.info(Header + "" + str(Intensity_sorted_mask) + "" + str(Ne_final_sorted_mask))  # Works well
            # logger.info(Intensity_sorted_mask, Ne_final_sorted_mask) #Works well

            FourPhotonAbsorption_params = curve_fit(FourPhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
            #ThreePhotonAbsorption_params = curve_fit(MultiPhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
            logger.info(Header+"** Ionization parameters for tau="+str(int(1e15*pulse_value))+" fs: "+str(FourPhotonAbsorption_params[0]))

            if(ShowPulseDurations):
                tau_text=r', $\tau=$'+str(int(pulse_value*1E15))+' fs'
            else:
                tau_text=""

            if(ShowSigmaRates==2):
                sigma_txt=r"$\sigma_4=$"+str('%.4g'%FourPhotonAbsorption_params[0][0])+tau_text
            elif(ShowSigmaRates==1):
                sigma_txt=r"$n=4$"
            else:
                sigma_txt=""

            logger.info(Header+"** Info: adding the fitted curve to the plot")
            Nexc_4PI = FourPhotonAbsorption(Intensity_sorted_mask, FourPhotonAbsorption_params[0])
            #print Nexc_2PI
            if(plotMe):
                plot_field, = ax1.plot(np.multiply(xscale,Intensity_sorted_mask), Nexc_4PI, Colors[pulse_value], linestyle=linestyle_tuple[LineType_MPIfitArray[4]], label=sigma_txt) #NOTE adds value of sigma_3 here
            else:
                plot_field = []
            del Nexc_4PI, FourPhotonAbsorption_params, Intensity_sorted_mask
            return plot_field
        #}}}

        def FitFivePhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, n_photons, FivePhotonIntensityRegime_min, FivePhotonIntensityRegime_max, plotMe=True): #TODO: pass function name in argument instead of copying this function... :-(
            #{{{
            logger.info(Header+"** Info: Filtering the fitting region...")
            Intensity_sorted_np = np.array(Intensity_sorted) #filtering works for NumPy array, not for lists.
            #condition = (Intensity_sorted_np > TwoPhotonIntensityRegime_min) and (Intensity_sorted_np < TwoPhotonIntensityRegime_max)
            Intensity_sorted_mask = Intensity_sorted_np[(Intensity_sorted_np > FivePhotonIntensityRegime_min) & (Intensity_sorted_np < FivePhotonIntensityRegime_max)]

            Ne_final_sorted_np  = np.array(Ne_final_sorted)
            Ne_final_sorted_mask  = Ne_final_sorted_np[(Intensity_sorted_np > FivePhotonIntensityRegime_min) & (Intensity_sorted_np < FivePhotonIntensityRegime_max)]
            # logger.info(Intensity_sorted_mask, Ne_final_sorted_mask) #Works well
            logger.info(Header + "" + str(Intensity_sorted_mask) + "" + str(Ne_final_sorted_mask))  # Works well

            FivePhotonAbsorption_params = curve_fit(FivePhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
            #ThreePhotonAbsorption_params = curve_fit(MultiPhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
            logger.info(Header+"** Ionization parameters for tau="+str(int(1e15*pulse_value))+" fs: "+str(FivePhotonAbsorption_params[0]))

            logger.info(Header+"** Info: adding the fitted curve to the plot")
            Nexc_5PI = FivePhotonAbsorption(Intensity_sorted_mask, FivePhotonAbsorption_params[0])
            #print Nexc_2PI

            if(ShowPulseDurations):
                tau_text=r', $\tau=$'+str(int(pulse_value*1E15))+' fs'
            else:
                tau_text=""

            if(ShowSigmaRates==2):
                sigma_txt=r"$\sigma_5=$"+str('%.4g'%FivePhotonAbsorption_params[0][0])+tau_text
            elif(ShowSigmaRates==1):
                sigma_txt=r"$n=5$"
            else:
                sigma_txt=""
            if(plotMe):
                plot_field, = ax1.plot(np.multiply(xscale,Intensity_sorted_mask), Nexc_5PI, Colors[pulse_value], linestyle=linestyle_tuple[LineType_MPIfitArray[5]], label=sigma_txt) #NOTE adds value of sigma_3 here
            else:
                plot_field = []
            del Nexc_5PI, FivePhotonAbsorption_params, Intensity_sorted_mask
            return plot_field
        #}}}

        def FitSixPhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, n_photons, SixPhotonIntensityRegime_min, SixPhotonIntensityRegime_max, plotMe=True): #TODO: pass function name in argument instead of copying this function... :-(
            #{{{
            logger.info(Header+"** Info: Filtering the fitting region...")
            Intensity_sorted_np = np.array(Intensity_sorted) #filtering works for NumPy array, not for lists.
            #condition = (Intensity_sorted_np > TwoPhotonIntensityRegime_min) and (Intensity_sorted_np < TwoPhotonIntensityRegime_max)
            Intensity_sorted_mask = Intensity_sorted_np[(Intensity_sorted_np > SixPhotonIntensityRegime_min) & (Intensity_sorted_np < SixPhotonIntensityRegime_max)]

            Ne_final_sorted_np  = np.array(Ne_final_sorted)
            Ne_final_sorted_mask  = Ne_final_sorted_np[(Intensity_sorted_np > SixPhotonIntensityRegime_min) & (Intensity_sorted_np < SixPhotonIntensityRegime_max)]
            # logger.info(Intensity_sorted_mask, Ne_final_sorted_mask) #Works well
            logger.info(Header + "" + str(Intensity_sorted_mask) + "" + str(Ne_final_sorted_mask))  # Works well

            SixPhotonAbsorption_params = curve_fit(SixPhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
            #ThreePhotonAbsorption_params = curve_fit(MultiPhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
            logger.info(Header+"** Ionization parameters for tau="+str(int(1e15*pulse_value))+" fs: "+str(SixPhotonAbsorption_params[0]))

            logger.info(Header+"** Info: adding the fitted curve to the plot")
            Nexc_6PI = SixPhotonAbsorption(Intensity_sorted_mask, SixPhotonAbsorption_params[0])
            #print Nexc_2PI

            if(ShowPulseDurations):
                tau_text=r', $\tau=$'+str(int(pulse_value*1E15))+' fs'
            else:
                tau_text=""

            if(ShowSigmaRates==2):
                sigma_txt=r"$\sigma_6=$"+str('%.4g'%SixPhotonAbsorption_params[0][0])+tau_text
            elif(ShowSigmaRates==1):
                sigma_txt=r"$n=6$"
            else:
                sigma_txt=""
            if(plotMe):
                plot_field, = ax1.plot(np.multiply(xscale,Intensity_sorted_mask), Nexc_6PI, Colors[pulse_value], linestyle=linestyle_tuple[LineType_MPIfitArray[6]], label=sigma_txt) #NOTE adds value of sigma_3 here
            else:
                plot_field = []
            del Nexc_6PI, SixPhotonAbsorption_params, Intensity_sorted_mask
            return plot_field
        #}}}

        def FitSevenPhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, n_photons, SevenPhotonIntensityRegime_min, SevenPhotonIntensityRegime_max, plotMe=True): #TODO: pass function name in argument instead of copying this function... :-(
            #{{{
            logger.info(Header+"** Info: Filtering the fitting region...")
            Intensity_sorted_np = np.array(Intensity_sorted) #filtering works for NumPy array, not for lists.
            #condition = (Intensity_sorted_np > TwoPhotonIntensityRegime_min) and (Intensity_sorted_np < TwoPhotonIntensityRegime_max)
            Intensity_sorted_mask = Intensity_sorted_np[(Intensity_sorted_np > SevenPhotonIntensityRegime_min) & (Intensity_sorted_np < SevenPhotonIntensityRegime_max)]

            Ne_final_sorted_np  = np.array(Ne_final_sorted)
            Ne_final_sorted_mask  = Ne_final_sorted_np[(Intensity_sorted_np > SevenPhotonIntensityRegime_min) & (Intensity_sorted_np < SevenPhotonIntensityRegime_max)]
            # logger.info(Intensity_sorted_mask, Ne_final_sorted_mask) #Works well
            logger.info(Header + "" + str(Intensity_sorted_mask) + "" + str(Ne_final_sorted_mask))  # Works well

            SevenPhotonAbsorption_params = curve_fit(SevenPhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
            #ThreePhotonAbsorption_params = curve_fit(MultiPhotonAbsorption, Intensity_sorted_mask, Ne_final_sorted_mask)
            logger.info(Header+"** Ionization parameters for tau="+str(int(1e15*pulse_value))+" fs: "+str(SevenPhotonAbsorption_params[0]))

            logger.info(Header+"** Info: adding the fitted curve to the plot")
            Nexc_7PI = SevenPhotonAbsorption(Intensity_sorted_mask, SevenPhotonAbsorption_params[0])
            #print Nexc_2PI

            if(ShowPulseDurations):
                tau_text=r', $\tau=$'+str(int(pulse_value*1E15))+' fs'
            else:
                tau_text=""

            if(ShowSigmaRates==2):
                sigma_txt=r"$\sigma_7=$"+str('%.4g'%SevenPhotonAbsorption_params[0][0])+tau_text
            elif(ShowSigmaRates==1):
                sigma_txt=r"$n=7$"
            else:
                sigma_txt=""

            if(plotMe):
                plot_field, = ax1.plot(np.multiply(xscale,Intensity_sorted_mask), Nexc_7PI, Colors[pulse_value], linestyle=linestyle_tuple[LineType_MPIfitArray[7]], label=sigma_txt) #NOTE adds value of sigma_3 here
            else:
                plot_field = []
            del Nexc_7PI, SevenPhotonAbsorption_params, Intensity_sorted_mask
            return plot_field
        #}}}

        #NOTE: the photon number does NOT change monotoneously with increasing intensity
        #TODO: one has to code a fitting method that sweeps intensity, but accepts ANY number of photons on each interval.
        #NOTE: Maybe use the computed error to accept or reject each fitting? One can try all numbers of photons on each interval.
        #TODO: this part is functional-dependent. It should be generalized (these numbers are fine for LDA, not for TB09).
        #OnePhotonIntensityRegime_min=0; OnePhotonIntensityRegime_max=0; 
        if(Functional[0]==libPlotting.IndexFunctional['LDA']):
            logger.info("Entered LDA functional fitting...")
            if(pulseWavelength==483e-9):
                if(pulse_value==30e-15):
                    SevenPhotonIntensityRegime_min = 0e0;  SevenPhotonIntensityRegime_max  = 0E0
                    SixPhotonIntensityRegime_min = 0e0;    SixPhotonIntensityRegime_max  = 0E0
                    FivePhotonIntensityRegime_min = 0e0;   FivePhotonIntensityRegime_max  = 0E0
                    FourPhotonIntensityRegime_min = 0e0;   FourPhotonIntensityRegime_max = 0e0 #W/m2
                    ThreePhotonIntensityRegime_min = 0e0;  ThreePhotonIntensityRegime_max = 0e0 #W/m2
                    TwoPhotonIntensityRegime_min = 0e0;    TwoPhotonIntensityRegime_max = 0E0 #W/m2
                    OnePhotonIntensityRegime_min = 1.22e9*1E4; OnePhotonIntensityRegime_max = 1.4e13*1E4 #W/m2
                elif(pulse_value==20e-15):
                    SevenPhotonIntensityRegime_min = 0e0;  SevenPhotonIntensityRegime_max  = 0E0
                    SixPhotonIntensityRegime_min = 0e0;    SixPhotonIntensityRegime_max  = 0E0
                    FivePhotonIntensityRegime_min = 0e0;   FivePhotonIntensityRegime_max  = 0E0
                    FourPhotonIntensityRegime_min = 0e0;   FourPhotonIntensityRegime_max = 0e0 #W/m2
                    ThreePhotonIntensityRegime_min = 0e0;  ThreePhotonIntensityRegime_max = 0e0 #W/m2
                    TwoPhotonIntensityRegime_min = 0e0;    TwoPhotonIntensityRegime_max = 0E0 #W/m2
                    OnePhotonIntensityRegime_min = 1.22e9*1E4; OnePhotonIntensityRegime_max = 1.4e13*1E4 #W/m2
                else: #10 fs
                    SevenPhotonIntensityRegime_min = 0e0;  SevenPhotonIntensityRegime_max  = 0E0
                    SixPhotonIntensityRegime_min = 0e0;    SixPhotonIntensityRegime_max  = 0E0
                    FivePhotonIntensityRegime_min = 0e0;   FivePhotonIntensityRegime_max  = 0E0
                    FourPhotonIntensityRegime_min = 0e0;   FourPhotonIntensityRegime_max = 0e0 #W/m2
                    ThreePhotonIntensityRegime_min = 0e0;  ThreePhotonIntensityRegime_max = 0e0 #W/m2
                    TwoPhotonIntensityRegime_min = 0e0;    TwoPhotonIntensityRegime_max = 0E0 #W/m2
                    OnePhotonIntensityRegime_min = 1.22e9*1E4; OnePhotonIntensityRegime_max = 1.4e13*1E4 #W/m2

                #ListOfInterestingGammas = [0.10, 1.0, 10.]

            elif(pulseWavelength==800e-9):
                if(pulse_value==30e-15):
                    SevenPhotonIntensityRegime_min = 0e0;  SevenPhotonIntensityRegime_max  = 0E0
                    SixPhotonIntensityRegime_min = 0e0;    SixPhotonIntensityRegime_max  = 0E0
                    FivePhotonIntensityRegime_min = 0e0;   FivePhotonIntensityRegime_max  = 0E0
                    FourPhotonIntensityRegime_min = 0e0;   FourPhotonIntensityRegime_max = 0e0 #W/m2
                    ThreePhotonIntensityRegime_min = 0e0;  ThreePhotonIntensityRegime_max = 0e0 #W/m2
                    TwoPhotonIntensityRegime_min = 1.22e9*1E4; TwoPhotonIntensityRegime_max = 1.4e11*1E4 #W/m2
                    OnePhotonIntensityRegime_min = 1.4e11*1E4; OnePhotonIntensityRegime_max = 1.4e13*1E4 #W/m2
                elif(pulse_value==20e-15):
                    SevenPhotonIntensityRegime_min = 0e0;  SevenPhotonIntensityRegime_max  = 0E0
                    SixPhotonIntensityRegime_min = 0e0;    SixPhotonIntensityRegime_max  = 0E0
                    FivePhotonIntensityRegime_min = 0e0;   FivePhotonIntensityRegime_max  = 0E0
                    FourPhotonIntensityRegime_min = 0e0;   FourPhotonIntensityRegime_max = 0e0 #W/m2
                    ThreePhotonIntensityRegime_min = 0e0;  ThreePhotonIntensityRegime_max = 0e0 #W/m2
                    TwoPhotonIntensityRegime_min = 1.22e9*1E4; TwoPhotonIntensityRegime_max = 1.4e11*1E4 #W/m2
                    OnePhotonIntensityRegime_min = 1.4e11*1E4; OnePhotonIntensityRegime_max = 1.4e13*1E4 #W/m2
                else: #10 fs
                    SevenPhotonIntensityRegime_min = 0e0;  SevenPhotonIntensityRegime_max  = 0E0
                    SixPhotonIntensityRegime_min = 0e0;    SixPhotonIntensityRegime_max  = 0E0
                    FivePhotonIntensityRegime_min = 0e0;   FivePhotonIntensityRegime_max  = 0E0
                    FourPhotonIntensityRegime_min = 0e0;   FourPhotonIntensityRegime_max = 0e0 #W/m2
                    ThreePhotonIntensityRegime_min = 0e0;  ThreePhotonIntensityRegime_max = 0e0 #W/m2
                    TwoPhotonIntensityRegime_min = 1.22e9*1E4; TwoPhotonIntensityRegime_max = 1.4e11*1E4 #W/m2
                    OnePhotonIntensityRegime_min = 1.4e11*1E4; OnePhotonIntensityRegime_max = 1.4e13*1E4 #W/m2

                #ListOfInterestingGammas = [0.10, 1.0, 10.]

            elif(pulseWavelength==1600e-9):
                if(pulse_value==30e-15):
                    SevenPhotonIntensityRegime_min= 0e0;    SevenPhotonIntensityRegime_max  = 0E0
                    SixPhotonIntensityRegime_min  = 0e0;    SixPhotonIntensityRegime_max  = 0E0
                    FivePhotonIntensityRegime_min = 0e0;    FivePhotonIntensityRegime_max  = 0E0
                    FourPhotonIntensityRegime_min = 1.1E14; FourPhotonIntensityRegime_max = 1.6E15 #W/m2
                    ThreePhotonIntensityRegime_min=3.39E15; ThreePhotonIntensityRegime_max  = 5.29E15
                    TwoPhotonIntensityRegime_min = 1.4E15;  TwoPhotonIntensityRegime_max = 3.2E16
                    OnePhotonIntensityRegime_min = 3.2E16;  OnePhotonIntensityRegime_max = 1.6E17
                elif(pulse_value==20e-15):
                    SevenPhotonIntensityRegime_min  = 0e0;  SevenPhotonIntensityRegime_max  = 0E0
                    SixPhotonIntensityRegime_min  = 0e0;    SixPhotonIntensityRegime_max  = 0E0
                    FivePhotonIntensityRegime_min = 0e0;    FivePhotonIntensityRegime_max  = 0E0
                    FourPhotonIntensityRegime_min = 1.1E14; FourPhotonIntensityRegime_max = 1.6E15 #W/m2
                    ThreePhotonIntensityRegime_min=3.39E15; ThreePhotonIntensityRegime_max  = 5.29E15
                    TwoPhotonIntensityRegime_min = 1.3E15;  TwoPhotonIntensityRegime_max = 3.2E16
                    OnePhotonIntensityRegime_min = 3.24E16; OnePhotonIntensityRegime_max = 1.33E17
                else:
                    SevenPhotonIntensityRegime_min  = 0e0;  SevenPhotonIntensityRegime_max  = 0E0
                    SixPhotonIntensityRegime_min  = 0e0;    SixPhotonIntensityRegime_max  = 0E0
                    FivePhotonIntensityRegime_min  = 0e0;   FivePhotonIntensityRegime_max  = 0E0
                    FourPhotonIntensityRegime_min = 1.1E14; FourPhotonIntensityRegime_max = 1.6E15 #W/m2
                    ThreePhotonIntensityRegime_min  = 0e0;  ThreePhotonIntensityRegime_max  = 0e0
                    TwoPhotonIntensityRegime_min = 1.4E15;  TwoPhotonIntensityRegime_max = 3.2E16
                    OnePhotonIntensityRegime_min = 3.2E16;  OnePhotonIntensityRegime_max = 1.6E17
            elif (pulseWavelength == 1030e-9):
                if (pulse_value == 30e-15):
                    SevenPhotonIntensityRegime_min = 0e0; SevenPhotonIntensityRegime_max = 0E0
                    SixPhotonIntensityRegime_min = 0e0; SixPhotonIntensityRegime_max = 0E0
                    FivePhotonIntensityRegime_min = 0e0; FivePhotonIntensityRegime_max = 0E0
                    FourPhotonIntensityRegime_min = 0e0; FourPhotonIntensityRegime_max = 0e0  # W/m2
                    ThreePhotonIntensityRegime_min = 0e0; ThreePhotonIntensityRegime_max = 0e0  # W/m2
                    TwoPhotonIntensityRegime_min = np.power(10e0,13.12293944); TwoPhotonIntensityRegime_max = np.power(10e0,16.02420264) #2
                    OnePhotonIntensityRegime_min = np.power(10e0,16.02420264); OnePhotonIntensityRegime_max = np.power(10e0,17.33862237) #1
                    # TwoPhotonIntensityRegime_min = np.power(10e0, 16.02420264); TwoPhotonIntensityRegime_max = np.power(10e0, 17.33862237) #1
                    # OnePhotonIntensityRegime_min = np.power(10e0,17.33862237); OnePhotonIntensityRegime_max = np.power(10e0,19.031422446) #0
                elif (pulse_value == 20e-15):
                    SevenPhotonIntensityRegime_min = 0e0;
                    SevenPhotonIntensityRegime_max = 0E0
                    SixPhotonIntensityRegime_min = 0e0;
                    SixPhotonIntensityRegime_max = 0E0
                    FivePhotonIntensityRegime_min = 0e0;
                    FivePhotonIntensityRegime_max = 0E0
                    FourPhotonIntensityRegime_min = 0e0; FourPhotonIntensityRegime_max = 0e0  # W/m2
                    ThreePhotonIntensityRegime_min = 0e0;
                    ThreePhotonIntensityRegime_max = 0e0  # W/m2
                    TwoPhotonIntensityRegime_min = np.power(10e0, 13.12293944);
                    TwoPhotonIntensityRegime_max = np.power(10e0, 16.02420264)  # 2
                    OnePhotonIntensityRegime_min = np.power(10e0, 16.02420264);
                    OnePhotonIntensityRegime_max = np.power(10e0, 17.33862237)  # 1
                    # TwoPhotonIntensityRegime_min = np.power(10e0, 16.02420264); TwoPhotonIntensityRegime_max = np.power(10e0, 17.33862237) #1
                elif(pulse_value==10e-15):
                    SevenPhotonIntensityRegime_min = 0e0;
                    SevenPhotonIntensityRegime_max = 0E0
                    SixPhotonIntensityRegime_min = 0e0;
                    SixPhotonIntensityRegime_max = 0E0
                    FivePhotonIntensityRegime_min = 0e0;
                    FivePhotonIntensityRegime_max = 0E0
                    FourPhotonIntensityRegime_min = 0e0; FourPhotonIntensityRegime_max = 0e0  # W/m2
                    ThreePhotonIntensityRegime_min = 0e0;
                    ThreePhotonIntensityRegime_max = 0e0  # W/m2
                    TwoPhotonIntensityRegime_min = np.power(10e0, 13.12293944);
                    TwoPhotonIntensityRegime_max = np.power(10e0, 16.02420264)  # 2
                    OnePhotonIntensityRegime_min = np.power(10e0, 16.02420264);
                    OnePhotonIntensityRegime_max = np.power(10e0, 17.33862237)  # 1
                else:
                    logger.warning("Branching to default case in LDA, 1030 nm")
                    SevenPhotonIntensityRegime_min = 0e0; SevenPhotonIntensityRegime_max = 0E0
                    SixPhotonIntensityRegime_min = 0e0;   SixPhotonIntensityRegime_max = 0E0
                    FivePhotonIntensityRegime_min = 0e0;  FivePhotonIntensityRegime_max = 0E0
                    FourPhotonIntensityRegime_min = 0e0;  FourPhotonIntensityRegime_max = 0e0  # W/m2
                    ThreePhotonIntensityRegime_min = 0e0; ThreePhotonIntensityRegime_max = 0e0  # W/m2
                    TwoPhotonIntensityRegime_min = np.power(10e0, 15.12293944);     TwoPhotonIntensityRegime_max = np.power(10e0, 16.59865948)  # 2
                    OnePhotonIntensityRegime_min = 0e0; OnePhotonIntensityRegime_max = 0e0;  # 1
                    # TwoPhotonIntensityRegime_min = np.power(10e0, 15.1229);     TwoPhotonIntensityRegime_max = np.power(10e0, 15.80504902)  # 2
                    # OnePhotonIntensityRegime_min = np.power(10e0, 15.80504902); OnePhotonIntensityRegime_max = np.power(10e0, 16.86325262)  # 1
                    # TwoPhotonIntensityRegime_min = np.power(10e0, 16.02420264); TwoPhotonIntensityRegime_max = np.power(10e0, 17.33862237) #1
                #TwoPhotonIntensityRegime_min = 3.4E15
                #TwoPhotonIntensityRegime_max = 4.7E16 #W/m2

                #OnePhotonIntensityRegime_min = 4.7E16
                #OnePhotonIntensityRegime_max = 1.3E17 #W/m2

                #ListOfInterestingGammas = [0.10, 1.0, 10.]

            elif(pulseWavelength==3200e-9): #added to Table
                if(pulse_value==20e-15):
                    SevenPhotonIntensityRegime_min = 0E0 #2.1E14
                    SevenPhotonIntensityRegime_max = 0E0 #1.3E15
                    SixPhotonIntensityRegime_min = 0E0;      SixPhotonIntensityRegime_max = 0E0
                    FivePhotonIntensityRegime_min = 2.1E14;  FivePhotonIntensityRegime_max = 9.6E14
                    FourPhotonIntensityRegime_min = 9.6E14;     FourPhotonIntensityRegime_max  = 2.6E15
                    ThreePhotonIntensityRegime_min = 2.6E15; ThreePhotonIntensityRegime_max = 6.4E15 #W/m2
                    TwoPhotonIntensityRegime_min = 6.4E15;  TwoPhotonIntensityRegime_max = 3.4E16 #W/m2
                    OnePhotonIntensityRegime_min = 3.41E16;  OnePhotonIntensityRegime_max = 1.0E13*1E4 #7.3E16 #W/m2
                elif(pulse_value==30e-15):
                    SevenPhotonIntensityRegime_min = 0E0 #2.1E14
                    SevenPhotonIntensityRegime_max = 0E0 #1.3E15
                    SixPhotonIntensityRegime_min = 0E0;      SixPhotonIntensityRegime_max = 0E0
                    FivePhotonIntensityRegime_min = 2.1E14;  FivePhotonIntensityRegime_max = 9.9E14 #1.3E15
                    FourPhotonIntensityRegime_min = 9.9E14;     FourPhotonIntensityRegime_max  = 2.58E15
                    ThreePhotonIntensityRegime_min = 2.58E15; ThreePhotonIntensityRegime_max = 5.31E15 #W/m2
                    TwoPhotonIntensityRegime_min = 5.31E15;  TwoPhotonIntensityRegime_max = 3.4E16 #W/m2
                    OnePhotonIntensityRegime_min = 3.41E16;  OnePhotonIntensityRegime_max = 1.0E13*1E4 #7.3E16 #W/m2
                else: #works for 10 fs and 10.xx fs then
                    SevenPhotonIntensityRegime_min= 0E0 #2.1E14
                    SevenPhotonIntensityRegime_max= 0E0 #1.3E15
                    SixPhotonIntensityRegime_min  = 0E0;     SixPhotonIntensityRegime_max  = 0E0
                    FivePhotonIntensityRegime_min = 2.1E14;  FivePhotonIntensityRegime_max = 1.07E15
                    FourPhotonIntensityRegime_min = 1.07E15;     FourPhotonIntensityRegime_max = 2.6E15
                    ThreePhotonIntensityRegime_min= 2.87E15; ThreePhotonIntensityRegime_max= 8.3E15 #W/m2
                    TwoPhotonIntensityRegime_min = 8.2E15;  TwoPhotonIntensityRegime_max = 3.4E16 #W/m2
                    OnePhotonIntensityRegime_min = 3.41E16;  OnePhotonIntensityRegime_max = 1.0E13*1E4 #7.3E16 #W/m2

                #ListOfInterestingGammas = [0.10, 1.0, 10.]
            else:
                SevenPhotonIntensityRegime_min = 0e0; SevenPhotonIntensityRegime_max  = 0E0
                SixPhotonIntensityRegime_min = 0e0;   SixPhotonIntensityRegime_max  = 0E0
                FivePhotonIntensityRegime_min = 0e0;  FivePhotonIntensityRegime_max  = 0E0
                FourPhotonIntensityRegime_min = 0e0;  FourPhotonIntensityRegime_max  = 1.1E14
                ThreePhotonIntensityRegime_min = 0e0; ThreePhotonIntensityRegime_max = 0e0 #W/m2
                TwoPhotonIntensityRegime_min = 0e0;   TwoPhotonIntensityRegime_max = 1.3e15 #W/m2
                OnePhotonIntensityRegime_min = 1e12;  OnePhotonIntensityRegime_max = 1.3e15 #W/m2

                #ListOfInterestingGammas = [0.10, 1.0, 10.]
        elif(Functional[0]==libPlotting.IndexFunctional['TB09']):
            logger.info("Entered TB09 functional fitting...")
            if(abs(pulseWavelength-407.58e-9) < 1E-9 ):
                logger.info("Entered TB09 functional fitting for 407 nm...")
                if(pulse_value==30e-15):
                    logger.info("Entered TB09 functional fitting for 407 nm, 30 fs...")
                    SevenPhotonIntensityRegime_min = 0e0;  SevenPhotonIntensityRegime_max  = 0E0
                    SixPhotonIntensityRegime_min = 0e0;    SixPhotonIntensityRegime_max  = 0E0
                    FivePhotonIntensityRegime_min = 0e0;   FivePhotonIntensityRegime_max  = 0E0
                    FourPhotonIntensityRegime_min = 0e0;   FourPhotonIntensityRegime_max = 0e0 #W/m2
                    ThreePhotonIntensityRegime_min = 0e0;  ThreePhotonIntensityRegime_max = 0e0 #W/m2
                    TwoPhotonIntensityRegime_min = 0e0;    TwoPhotonIntensityRegime_max = 0E0 #W/m2
                    OnePhotonIntensityRegime_min = 1.22e9*1E4; OnePhotonIntensityRegime_max = 1.4e13*1E4 #W/m2
                elif(pulse_value==20e-15):
                    logger.info("Entered TB09 functional fitting for 407 nm, 20 fs...")
                    SevenPhotonIntensityRegime_min = 0e0;  SevenPhotonIntensityRegime_max  = 0E0
                    SixPhotonIntensityRegime_min = 0e0;    SixPhotonIntensityRegime_max  = 0E0
                    FivePhotonIntensityRegime_min = 0e0;   FivePhotonIntensityRegime_max  = 0E0
                    FourPhotonIntensityRegime_min = 0e0;   FourPhotonIntensityRegime_max = 0e0 #W/m2
                    ThreePhotonIntensityRegime_min = 0e0;  ThreePhotonIntensityRegime_max = 0e0 #W/m2
                    TwoPhotonIntensityRegime_min = 0e0;    TwoPhotonIntensityRegime_max = 0E0 #W/m2
                    OnePhotonIntensityRegime_min = 1.22e9*1E4; OnePhotonIntensityRegime_max = 1.4e13*1E4 #W/m2
                else: #10 fs
                    logger.info("Entered TB09 functional fitting for 407 nm, any fs...")
                    SevenPhotonIntensityRegime_min = 0e0;  SevenPhotonIntensityRegime_max  = 0E0
                    SixPhotonIntensityRegime_min = 0e0;    SixPhotonIntensityRegime_max  = 0E0
                    FivePhotonIntensityRegime_min = 0e0;   FivePhotonIntensityRegime_max  = 0E0
                    FourPhotonIntensityRegime_min = 0e0;   FourPhotonIntensityRegime_max = 0e0 #W/m2
                    ThreePhotonIntensityRegime_min = 0e0;  ThreePhotonIntensityRegime_max = 0e0 #W/m2
                    TwoPhotonIntensityRegime_min = 0e0;    TwoPhotonIntensityRegime_max = 0E0 #W/m2
                    OnePhotonIntensityRegime_min = 1.22e9*1E4; OnePhotonIntensityRegime_max = 1.4e13*1E4 #W/m2

                #ListOfInterestingGammas = [0.10, 1.0, 10.]

            elif(abs(pulseWavelength)-800e-9 < 1E-9):
                logger.info("Entered TB09 functional fitting for 800 nm...")
                if(pulse_value==30e-15):
                    logger.info("Entered TB09 functional fitting for 800 nm, 30 fs ...")
                    SevenPhotonIntensityRegime_min = 0e0;  SevenPhotonIntensityRegime_max  = 0E0
                    SixPhotonIntensityRegime_min = 0e0;    SixPhotonIntensityRegime_max  = 0E0
                    FivePhotonIntensityRegime_min = 0e0;   FivePhotonIntensityRegime_max  = 0E0
                    FourPhotonIntensityRegime_min = 0e0;   FourPhotonIntensityRegime_max = 0e0 #W/m2
                    ThreePhotonIntensityRegime_min = 0e0;  ThreePhotonIntensityRegime_max = 0e0 #W/m2
                    TwoPhotonIntensityRegime_min = 1.22e9*1E4; TwoPhotonIntensityRegime_max = 1.4e11*1E4 #W/m2
                    OnePhotonIntensityRegime_min = 1.4e11*1E4; OnePhotonIntensityRegime_max = 1.4e13*1E4 #W/m2
                elif(pulse_value==20e-15):
                    logger.info("Entered TB09 functional fitting for 800 nm, 20 fs ...")
                    SevenPhotonIntensityRegime_min = 0E0*1E4;  SevenPhotonIntensityRegime_max        = 0E0*1E4
                    SixPhotonIntensityRegime_min   = 0E0*1E4;    SixPhotonIntensityRegime_max        = 0E0*1E4
                    FivePhotonIntensityRegime_min  = 0E0*1E4;   FivePhotonIntensityRegime_max        = 0E0*1E4
                    FourPhotonIntensityRegime_min  = 0E0*1E4;   FourPhotonIntensityRegime_max        = 0E0*1E4 #W/m2
                    ThreePhotonIntensityRegime_min = 0E0*1E4;  ThreePhotonIntensityRegime_max = 0E0*1E4 #W/m2
                    TwoPhotonIntensityRegime_min   = 2.87e9*1E4; TwoPhotonIntensityRegime_max    = 1.21E12*1E4 #W/m2
                    OnePhotonIntensityRegime_min   = 1.60e12*1E4; OnePhotonIntensityRegime_max   = 1.34E13*1E4 #W/m2
                    #ThreePhotonIntensityRegime_min = 1.23e11*1E4;  ThreePhotonIntensityRegime_max = 1.6512*1E4 #W/m2
                    #TwoPhotonIntensityRegime_min = 2.99e9*1E4; TwoPhotonIntensityRegime_max = 1.23e11*1E4 #W/m2
                    #OnePhotonIntensityRegime_min = 1.23e11*1E4; OnePhotonIntensityRegime_max = 1.45E13*1E4 #W/m2
                else: #10 fs
                    logger.info("Entered TB09 functional fitting for 800 nm, any fs ...")
                    SevenPhotonIntensityRegime_min = 0e0;  SevenPhotonIntensityRegime_max  = 0E0
                    SixPhotonIntensityRegime_min = 0e0;    SixPhotonIntensityRegime_max  = 0E0
                    FivePhotonIntensityRegime_min = 0e0;   FivePhotonIntensityRegime_max  = 0E0
                    FourPhotonIntensityRegime_min = 0e0;   FourPhotonIntensityRegime_max = 0e0 #W/m2
                    ThreePhotonIntensityRegime_min = 0e0;  ThreePhotonIntensityRegime_max = 0e0 #W/m2
                    TwoPhotonIntensityRegime_min = 1.22e9*1E4; TwoPhotonIntensityRegime_max = 1.4e11*1E4 #W/m2
                    OnePhotonIntensityRegime_min = 1.4e11*1E4; OnePhotonIntensityRegime_max = 1.4e13*1E4 #W/m2
                #ListOfInterestingGammas = [0.10, 1.0, 10.]

            elif(abs(pulseWavelength-1600e-9)<1E-9):
                if(abs(pulse_value-30e-15))<1E-15:
                    SevenPhotonIntensityRegime_min= 0e0;    SevenPhotonIntensityRegime_max  = 0E0
                    SixPhotonIntensityRegime_min  = 0e0;    SixPhotonIntensityRegime_max  = 0E0
                    FivePhotonIntensityRegime_min = 0e0;    FivePhotonIntensityRegime_max  = 0E0
                    FourPhotonIntensityRegime_min = 1.1E14; FourPhotonIntensityRegime_max = 1.6E15 #W/m2
                    ThreePhotonIntensityRegime_min=3.39E15; ThreePhotonIntensityRegime_max  = 5.29E15
                    TwoPhotonIntensityRegime_min = 1.4E15;  TwoPhotonIntensityRegime_max = 3.2E16
                    OnePhotonIntensityRegime_min = 3.2E16;  OnePhotonIntensityRegime_max = 1.6E17
                elif(pulse_value==20e-15):
                    SevenPhotonIntensityRegime_min  = 0e0;  SevenPhotonIntensityRegime_max  = 0E0
                    SixPhotonIntensityRegime_min  = 0e0;    SixPhotonIntensityRegime_max  = 0E0
                    FivePhotonIntensityRegime_min = 0e0;    FivePhotonIntensityRegime_max  = 0E0
                    FourPhotonIntensityRegime_min = 1.1E14; FourPhotonIntensityRegime_max = 1.6E15 #W/m2
                    ThreePhotonIntensityRegime_min=3.39E15; ThreePhotonIntensityRegime_max  = 5.29E15
                    TwoPhotonIntensityRegime_min = 1.3E15;  TwoPhotonIntensityRegime_max = 3.2E16
                    OnePhotonIntensityRegime_min = 3.24E16; OnePhotonIntensityRegime_max = 1.33E17
                else:
                    SevenPhotonIntensityRegime_min  = 0e0;  SevenPhotonIntensityRegime_max  = 0E0
                    SixPhotonIntensityRegime_min  = 0e0;    SixPhotonIntensityRegime_max  = 0E0
                    FivePhotonIntensityRegime_min  = 0e0;   FivePhotonIntensityRegime_max  = 0E0
                    FourPhotonIntensityRegime_min = 1.1E14; FourPhotonIntensityRegime_max = 1.6E15 #W/m2
                    ThreePhotonIntensityRegime_min  = 0e0;  ThreePhotonIntensityRegime_max  = 0e0
                    TwoPhotonIntensityRegime_min = 1.4E15;  TwoPhotonIntensityRegime_max = 3.2E16
                    OnePhotonIntensityRegime_min = 3.2E16;  OnePhotonIntensityRegime_max = 1.6E17

                #TwoPhotonIntensityRegime_min = 3.4E15
                #TwoPhotonIntensityRegime_max = 4.7E16 #W/m2

                #OnePhotonIntensityRegime_min = 4.7E16
                #OnePhotonIntensityRegime_max = 1.3E17 #W/m2

                #ListOfInterestingGammas = [0.10, 1.0, 10.]

            elif(pulseWavelength==3200e-9): #added to Table
                if(pulse_value==20e-15):
                    SevenPhotonIntensityRegime_min = 0E0 #2.1E14
                    SevenPhotonIntensityRegime_max = 0E0 #1.3E15
                    SixPhotonIntensityRegime_min = 0E0;      SixPhotonIntensityRegime_max = 0E0
                    FivePhotonIntensityRegime_min = 2.1E14;  FivePhotonIntensityRegime_max = 9.6E14
                    FourPhotonIntensityRegime_min = 9.6E14;     FourPhotonIntensityRegime_max  = 2.6E15
                    ThreePhotonIntensityRegime_min = 2.6E15; ThreePhotonIntensityRegime_max = 6.4E15 #W/m2
                    TwoPhotonIntensityRegime_min = 6.4E15;  TwoPhotonIntensityRegime_max = 3.4E16 #W/m2
                    OnePhotonIntensityRegime_min = 3.41E16;  OnePhotonIntensityRegime_max = 1.0E13*1E4 #7.3E16 #W/m2
                elif(pulse_value==30e-15):
                    SevenPhotonIntensityRegime_min = 0E0 #2.1E14
                    SevenPhotonIntensityRegime_max = 0E0 #1.3E15
                    SixPhotonIntensityRegime_min = 0E0;      SixPhotonIntensityRegime_max = 0E0
                    FivePhotonIntensityRegime_min = 2.1E14;  FivePhotonIntensityRegime_max = 9.9E14 #1.3E15
                    FourPhotonIntensityRegime_min = 9.9E14;     FourPhotonIntensityRegime_max  = 2.58E15
                    ThreePhotonIntensityRegime_min = 2.58E15; ThreePhotonIntensityRegime_max = 5.31E15 #W/m2
                    TwoPhotonIntensityRegime_min = 5.31E15;  TwoPhotonIntensityRegime_max = 3.4E16 #W/m2
                    OnePhotonIntensityRegime_min = 3.41E16;  OnePhotonIntensityRegime_max = 1.0E13*1E4 #7.3E16 #W/m2
                else: #works for 10 fs and 10.xx fs then
                    SevenPhotonIntensityRegime_min= 0E0 #2.1E14
                    SevenPhotonIntensityRegime_max= 0E0 #1.3E15
                    SixPhotonIntensityRegime_min  = 0E0;     SixPhotonIntensityRegime_max  = 0E0
                    FivePhotonIntensityRegime_min = 2.1E14;  FivePhotonIntensityRegime_max = 1.07E15
                    FourPhotonIntensityRegime_min = 1.07E15;     FourPhotonIntensityRegime_max = 2.6E15
                    ThreePhotonIntensityRegime_min= 2.87E15; ThreePhotonIntensityRegime_max= 8.3E15 #W/m2
                    TwoPhotonIntensityRegime_min = 8.2E15;  TwoPhotonIntensityRegime_max = 3.4E16 #W/m2
                    OnePhotonIntensityRegime_min = 3.41E16;  OnePhotonIntensityRegime_max = 1.0E13*1E4 #7.3E16 #W/m2
            elif (pulseWavelength==1030E-9):
                logger.warning("Branching to 1030 nm case in TB09")
                SevenPhotonIntensityRegime_min = 0e0; SevenPhotonIntensityRegime_max = 0E0
                SixPhotonIntensityRegime_min = 0e0;   SixPhotonIntensityRegime_max = 0E0
                FivePhotonIntensityRegime_min = 0e0;  FivePhotonIntensityRegime_max = 0E0
                FourPhotonIntensityRegime_min = 0e0;  FourPhotonIntensityRegime_max = 0e0  # W/m2
                ThreePhotonIntensityRegime_min = 0e0; ThreePhotonIntensityRegime_max = 0e0  # W/m2
                TwoPhotonIntensityRegime_min = np.power(10e0, 15.12293944);     TwoPhotonIntensityRegime_max = np.power(10e0, 16.59865948)  # 2
                OnePhotonIntensityRegime_min = 0e0; OnePhotonIntensityRegime_max = 0e0;  # 1
                #ListOfInterestingGammas = [0.10, 1.0, 10.]
            else:
                SevenPhotonIntensityRegime_min = 0e0; SevenPhotonIntensityRegime_max  = 0E0
                SixPhotonIntensityRegime_min = 0e0;   SixPhotonIntensityRegime_max  = 0E0
                FivePhotonIntensityRegime_min = 0e0;  FivePhotonIntensityRegime_max  = 0E0
                FourPhotonIntensityRegime_min = 0e0;  FourPhotonIntensityRegime_max  = 1.1E14
                ThreePhotonIntensityRegime_min = 0e0; ThreePhotonIntensityRegime_max = 0e0 #W/m2
                TwoPhotonIntensityRegime_min = 0e0;   TwoPhotonIntensityRegime_max = 1.3e15 #W/m2
                OnePhotonIntensityRegime_min = 1e12;  OnePhotonIntensityRegime_max = 1.3e15 #W/m2
        else: #any other functional
            logger.info("Entered unknown functional analysis...")
            logger.debug("Functional contains: "+str(Functional))
            SevenPhotonIntensityRegime_min = 1e12; SevenPhotonIntensityRegime_max  = 1.3e15
            SixPhotonIntensityRegime_min   = 1e12;   SixPhotonIntensityRegime_max  = 1.3e15
            FivePhotonIntensityRegime_min  = 1e12;  FivePhotonIntensityRegime_max  = 1.3e15
            FourPhotonIntensityRegime_min  = 1e12;  FourPhotonIntensityRegime_max  = 1.3e15
            ThreePhotonIntensityRegime_min = 1e12; ThreePhotonIntensityRegime_max  = 1.3e15 #W/m2
            TwoPhotonIntensityRegime_min   = 1e12;   TwoPhotonIntensityRegime_max  = 1.3e15 #W/m2
            OnePhotonIntensityRegime_min   = 1e12;  OnePhotonIntensityRegime_max   = 1.3e15 #W/m2

        ListOfInterestingGammas = [0.10, 1.0, 10.]
        omega = 2.*np.pi*c/wavelengths1[0]
        NumberOfPhotons = np.divide(Egap_SI/hbar, omega)
        logger.info(NumberOfPhotons)
        NumberOfPhotons_round = np.ceil(NumberOfPhotons)
        logger.info(Header+"** INFO: NUMBER OF PHOTONS = "+str(NumberOfPhotons)+", rounded to: "+str(NumberOfPhotons_round))

        # Here we attempt to fit the excitation coefficient in the multiphotonic regime.
        # Hence, it depends on the number of photons required for overcoming the direct band gap energy.
        # NOTE: we employ only ONE band gap energy here, despite there is a full band structure to consider.
        # Note that tunneling effect only depends on the intensity. Therefore it could look like
        # 1-photon-absorption.
        PlotTheFittedLines=(ShowSigmaRates>0) #plots the fitting lines if ShowSigmaRates!=0
        try:
            #NOTE: assumption here is the MAIN BAND GAP is employed for the transition. But it may NOT be the case.
            #plotMe = (pulse_value == 20E-15)
            #if(NumberOfPhotons_round >= 4 and NumberOfPhotons_round < 5):
            logger.info(Header+"** Info: fitting with 7-photon absortion...")
            plot_field = FitSevenPhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, 7, SevenPhotonIntensityRegime_min, SevenPhotonIntensityRegime_max, PlotTheFittedLines)
            if(ShowSigmaRates>0):
                plot1 += [plot_field]
        except (NameError, ValueError):
            logger.warning(Header+"** Fit using 7-photon absorption failed. ")

        try:
            #NOTE: assumption here is the MAIN BAND GAP is employed for the transition. But it may NOT be the case.
            #plotMe = (pulse_value == 20E-15)
            #if(NumberOfPhotons_round >= 4 and NumberOfPhotons_round < 5):
            logger.info(Header+"** Info: fitting with 6-photon absortion...")
            plot_field = FitSixPhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, 6, SixPhotonIntensityRegime_min, SixPhotonIntensityRegime_max, PlotTheFittedLines)
            if(ShowSigmaRates>0):
                plot1 += [plot_field]
        except (NameError, ValueError):
            logger.warning(Header+"** Fit using 6-photon absorption failed. ")

        try:
            #NOTE: assumption here is the MAIN BAND GAP is employed for the transition. But it may NOT be the case.
            #plotMe = (pulse_value == 20E-15)
            #if(NumberOfPhotons_round >= 4 and NumberOfPhotons_round < 5):
            logger.info(Header+"** Info: fitting with 5-photon absortion...")
            plot_field = FitFivePhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, 5, FivePhotonIntensityRegime_min, FivePhotonIntensityRegime_max, PlotTheFittedLines)
            if(ShowSigmaRates>0):
                plot1 += [plot_field]
        except (NameError, ValueError):
            logger.warning(Header+"** Fit using 5-photon absorption failed. ")

        try:
            #NOTE: assumption here is the MAIN BAND GAP is employed for the transition. But it may NOT be the case.
            #plotMe = (pulse_value == 20E-15)
            #if(NumberOfPhotons_round >= 4 and NumberOfPhotons_round < 5):
            logger.info(Header+"** Info: fitting with 4-photon absortion...")
            plot_field = FitFourPhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, 4, FourPhotonIntensityRegime_min, FourPhotonIntensityRegime_max, PlotTheFittedLines)
            if(ShowSigmaRates>0):
                plot1 += [plot_field]
        except (NameError, ValueError):
            logger.warning(Header+"** Fit using 4-photon absorption failed. ")
        try:
            #elif(NumberOfPhotons_round >= 3 and NumberOfPhotons_round < 4):
            print(Header+"** Info: fitting with 3-photon absortion...")
            plot_field = FitThreePhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, 3, ThreePhotonIntensityRegime_min, ThreePhotonIntensityRegime_max, PlotTheFittedLines)
            if(ShowSigmaRates>0):
                plot1 += [plot_field]
        except:
            logger.warning(Header+"** Fit using 3-photon absorption failed. ")
        try:
            #elif(NumberOfPhotons_round >= 2 and NumberOfPhotons_round < 3):
            logger.info(Header+"** Info: fitting with 2-photon absortion...")
            plot_field = FitTwoPhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, 2, TwoPhotonIntensityRegime_min, TwoPhotonIntensityRegime_max, PlotTheFittedLines)
            if(ShowSigmaRates>0):
                plot1 += [plot_field]
        except (NameError, ValueError):
            logger.warning(Header+"** Fit using 2-photon absorption failed. ")
        try:
        #elif (NumberOfPhotons_round >= 1 and NumberOfPhotons_round < 2):
            logger.info(Header+"** Info: fitting with 1-photon absortion...")
            plot_field = FitOnePhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, 1, OnePhotonIntensityRegime_min, OnePhotonIntensityRegime_max, PlotTheFittedLines)
            if(ShowSigmaRates>0):
                plot1 += [plot_field]
        except (NameError, ValueError):
            logger.warning(Header+"** Fit using 1-photon absorption failed. ")
            #else:
            #print Header+"** Warning: could not fit anything, surely because you are using the 2-band vision! :) "

        #except:
        #print Header+"** Fit using N-photon absorption failed. "
        # NOTE: this is NOT 3-photon absorption regime, this is NON-PERTURBATIVE regime. Transition to tunneling.

        logger.info(Header+"** Info: fitting the non-perturbative regime...")
        TunnelingIntensityRegime_min = 2E15 #TwoPhotonIntensityRegime_max #1e12
        TunnelingIntensityRegime_max = 2E17  #1.3e15 #W/m2
        #try:
        #print Header+"** Info: fitting the tunneling rate..."
        #plot_field = FitOnePhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, 1, TunnelingIntensityRegime_min, TunnelingIntensityRegime_max)
        #plot1 += [plot_field]
        #except:
        #print Header+"** Fit using tunneling rate failed. "
        #try:
        #plot_field = FitTwoPhotonIonizationRateFromNeIntensity(Intensity_sorted, Ne_final_sorted, 2, ThreePhotonIntensityRegime_min, ThreePhotonIntensityRegime_max)
        #plot1 += [plot_field]
        #except:
        #print Header+"** Warning: 3-ph fitting of the data failed for tau="+str(pulse_value*1E15)+" fs."
        for gammaValues in ListOfInterestingGammas:
            ax1.axvline(x=np.multiply(xscale, libKeldysh.IntensityAtGamma(gammaValues, Egap_SI, meff_Waroquiers_LDA, pulseWavelength)), color='k', linestyle='--', linewidth=0.5)
        ax1.axhline(y=  Density, color='k', linestyle='--', linewidth=0.25)
        #ax1.axhline(y=2*Density, color='k', linestyle='--', linewidth=0.25)
        #ax1.axhline(y=3*Density, color='k', linestyle='--', linewidth=0.25)
        ax1.axhline(y=4*Density, color='k', linestyle='--', linewidth=0.25)
    #}}}
    ax1.set_xlabel(r"Intensity (W/cm$^{2}$)")
    #ax1.set_ylabel(r"Final exc. density (m$^{-3}$)")
    ax1.set_ylabel(r"n$_{\mathrm{exc}}$ (m$^{-3}$)")
    #plt.xlabel(r'Peak intensity (W/m$^2$)')
    #plt.ylabel(r'Final exc. density (m$^{-3}$)')
    #ax1.set_xscale('log')
    #ax1.set_yscale('log')
    #ax1.yaxis.set_minor_locator(mp.ticker.LogLocator(base=10.0, subs='auto', numticks=None))
    # TODO: Solution to use minor ticks is given here: https://stackoverflow.com/questions/44078409/matplotlib-semi-log-plot-minor-tick-marks-are-gone-when-range-is-large
    labels1 = [l.get_label() for l in plot1]
    ax1.legend(plot1, labels1, loc='best', fancybox=True, framealpha=1, fontsize=OneSizeFitsAll, ncol=ncols)

    #plt.tight_layout()
    #print Header+"Adding the adjusted coefficients to check our fitting law."

    ax2 = plt.subplot(plotB)
    if(MoveIon):
        MoveIonTitle=", w/ phonons."
        MoveIonFile="MoveIon"
    else:
        MoveIonTitle=""
        MoveIonFile=""
    if(SidePanels):
        ax2.set_title(r'$\lambda = $'+str(int(pulseWavelength*1E9))+r' nm, $n=\left\lceil\frac{E_g^\Gamma}{\hbar \omega}\right\rceil=$'+str(int(NumberOfPhotons_round))+MoveIonTitle)
    else:
        ax1.set_title(r'$\lambda = $'+str(int(pulseWavelength*1E9))+r' nm, $n=\left\lceil\frac{E_g^\Gamma}{\hbar \omega}\right\rceil=$'+str(int(NumberOfPhotons_round))+MoveIonTitle)
    logger.info(Header+"Generating comparison with Keldysh-Gruzdev model...") ## Looking forward to plot w_PI[intensity] from Keldysh_Gruzdev theory.
    order = 50
    try: 
        Efieldmin = np.min(Efield1_uniq)
        Efieldmax = np.max(Efield1_uniq) #5E10
    except:
        logger.error("Could not prepare the run for wavelength = "+str(pulseWavelength)+" m, Efield1_uniq = "+str(Efield1_uniq)+". It can be linked to hard filtering of pulse duration in program.")
        Efieldmin = 1E9
        Efieldmax = 1E11
    Efield1_Keldysh_log = np.linspace(np.log10(Efieldmin), np.log10(Efieldmax), 2000)
    Efield1_Keldysh     = np.power(10., Efield1_Keldysh_log)
    logger.debug(Efield1_Keldysh)

    for gammaValues in ListOfInterestingGammas:
        ax2.axvline(x=xscale*libKeldysh.IntensityAtGamma(gammaValues, Egap_SI, meff_Waroquiers_LDA, pulseWavelength), color='k', linestyle='--', linewidth=0.5)
    #  plt.ylim((1E16,4.*5e28))
    #ax1.grid()
    #if(PlotWpiAxis==True): #{{{
    plot2=[]


    if(EnableInsetMagnification):
        if(SidePanels):
            ax22 = inset_axes(ax2, width="100%", height="100%", borderpad=1.05, bbox_to_anchor=(0.65, 0.5, 1-0.65, 0.27), bbox_transform=ax2.transAxes, loc="center right")

            ax22.tick_params(top=False, left=False, bottom=True, labelleft=False, labeltop=False, labelbottom=True)
        else: #Panels one below the another
            ax22 = inset_axes(ax2, width="100%", height="100%", borderpad=1.05, bbox_to_anchor=(0.05, 0.6, 0.20, 1-0.6), bbox_transform=ax2.transAxes, loc="center right") #, zorder=6) #NOTE: only from Python3
            ax22.tick_params(top=False, left=False, bottom=True, labelleft=False, labeltop=False, labelbottom=True)

    for pulse_value in pulses_uniq: #{{{
        # select the corresponding wavelength in results
        PlottingSinglePulses_pulses = libDatabase.FilterDatabase(PlottingSinglePulses, pulse_value, libPlotting.IndexTable['pulses1'])
        PlottingSinglePulses_pulses = libDatabase.FilterDatabase(PlottingSinglePulses_pulses, pulse_value, libPlotting.IndexTable['pulses2']) #keeps only pulses of the same duration.

        runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom, LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingTable)

        if(UseNexcAfterPulse):
            del Ne_final
            Ne_final = Ne_afterpulse #NOTE: dirty and valid only until next call of ExtractData

        if(ShowPulseDurations):
            tau_text=r'$\tau=$'+str(int(pulse_value*1E15))+' fs'
        else:
            tau_text=""

        lists = sorted(zip(*[libKeldysh.FieldToIntensity(Efield1), Ne_final, pulses1]))
        Intensity_sorted, Ne_final_sorted, pulses1_sorted = list(zip(*lists))


        plot21, = ax2.plot(np.multiply(xscale, Intensity_sorted), np.divide(Ne_final_sorted, np.multiply(1.,pulses1_sorted)), Colors[pulse_value]+Symbols[pulse_value], linewidth=LineWidth[pulse_value], label=r'$w_{\mathrm{PI}}^{\mathrm{TDDFT}}$'+tau_text)
        plot2 += [plot21]
        if(EnableInsetMagnification):
            plot212, = ax22.plot(np.multiply(xscale, Intensity_sorted), np.divide(Ne_final_sorted, np.multiply(1.,pulses1_sorted)), Colors[pulse_value]+Symbols[pulse_value], linewidth=LineWidth[pulse_value], label=r'$w_{\mathrm{PI}}^{\mathrm{TDDFT}}$'+tau_text)

        #plot21, = ax2.plot(np.multiply(xscale, Intensity_sorted), Saturation, Colors[pulse_value]+Symbols[pulse_value], label=r'$\tilde{w}_{\mathrm{PI}}^{\mathrm{sat.}}$, $\tau=$'+str(int(pulse_value*1E15))+' fs')
        #plot2 += [plot21]

    #}}}

    if(pulseWavelength==800e-9):
        EffectiveMassList=[meff_Waroquiers_LDA] #[0.05]
        #EffectiveMassList=[meff_Waroquiers_LDA, 1.] #, meff_Waroquiers_LDA]
        #EnableInsetMagnification=False
    elif (pulseWavelength == 1030e-9):
        EffectiveMassList = [meff_Waroquiers_LDA]
    elif(pulseWavelength==1600e-9):
        EffectiveMassList=[meff_Waroquiers_LDA]
        #EnableInsetMagnification=False
    elif(pulseWavelength==3200e-9):
        EffectiveMassList=[meff_Waroquiers_LDA] #[0.025]
    #, meff_Waroquiers_LDA]
    #EnableInsetMagnification=True
    else:
        EffectiveMassList=[0.01, 0.025, 0.05, 0.1, meff_Waroquiers_LDA]

    for meff_Keldysh in EffectiveMassList: #{{{
        PeakWpi_Gruzdev = libKeldysh.GenerateKeldyshDatabase(Egap_SI, meff_Keldysh, pulseWavelength, Efield1_Keldysh*np.sqrt(OpticalIndex), order, 1.) #NOTE: optical index should be treated as such, not through adiabadicity coefficient, or?
        Saturation = (4.*Ntotal - (PeakWpi_Gruzdev * 20E-15 )) / (4.*Ntotal) #saturation term for a given pulse duration. of course, it depends on pulse duration ! For now, let's apply it only for 20 fs. Otherwise, it would be too heavy.
        PeakWpi_Gruzdev_sat = PeakWpi_Gruzdev * Saturation

        Acoeff=1
        PeakWpi_Atomic_Keldysh = libKeldysh.IonizationRateAtoms(pulseWavelength, meff_Keldysh, Egap_SI, Efield1_Keldysh, order, 1., Ntotal) #Optical index in gas in 1.
        logger.info(Header+"Info: Atomic Keldysh: "+str(np.min(PeakWpi_Atomic_Keldysh))+", "+str(np.max(PeakWpi_Atomic_Keldysh)))
        #PeakWpi_Gulley = GenerateKeldyshGulleyDatabase(Egap_SI, meff_Keldysh, pulseWavelength, Efield1_Keldysh, order, 1.)
        gamma = libKeldysh.gammaKeldysh(Egap_SI, meff_Keldysh, Efield1_Keldysh, pulseWavelength, 1.) #for gamma, ONE MUST NOT USE OPTICAL INDEX (because gamma is oscillation time divided by tunneling time, hence optical index must vanish)
        #ax2x2 = ax2.twiny()
        # We need to plot values of intensities for which gamma=0.1, 1, 3, 5, 10.
        # d e f IntensityAtGamma(gamma, Egap_SI, meff, wavelength):
        #omega = 2.*np.pi * c / wavelength
        #Ipeak = meff*m_e*Egap_SI*omega**2*c*epsilon_0 / (2. * e**2 * gamma ** 2)
        #return Ipeak

        ## GRUZDEV
        ## NOTE: Gruzdev treak consists in plotting (external_field, w_pi[internal_field]).
        plot22, = ax2.plot(np.multiply(xscale, libKeldysh.FieldToIntensity(Efield1_Keldysh)), PeakWpi_Gruzdev, '-k', label=r'$w_{\mathrm{PI}}^{\mathrm{KG}}$') #, $m_{\mathrm{eff}}=$'+str(np.round(meff_Keldysh,3))+'.')
        plot2    += [plot22]

        if(EnableInsetMagnification):
            plot222, = ax22.plot(np.multiply(xscale, libKeldysh.FieldToIntensity(Efield1_Keldysh)), PeakWpi_Gruzdev, '-k', label=r'$w_{\mathrm{PI}}^{\mathrm{KG}}$')

        if(PlotSaturation):
            plot22_sat, = ax2.plot(np.multiply(xscale, libKeldysh.FieldToIntensity(Efield1_Keldysh)), PeakWpi_Gruzdev_sat, '--.k', label=r'$w_{\mathrm{PI}}^{\mathrm{KG}} \times \frac{Z \times n_0-n_{\mathrm{exc}}}{Z \times n_0}$') #, $m_{\mathrm{eff}}=$'+str(np.round(meff_Keldysh,3))+'.')
            plot2    += [plot22_sat]

            if(EnableInsetMagnification):
                plot222_sat, = ax22.plot(np.multiply(xscale, libKeldysh.FieldToIntensity(Efield1_Keldysh)), PeakWpi_Gruzdev_sat, '--.k', label=r'$w_{\mathrm{PI}}^{\mathrm{KG}} \times \frac{Z \times n_0-n_{\mathrm{exc}}}{Z \times n_0}$')

        # ATOMIC KELDYSH
        logger.info(Header+"Preparing atomic Keldysh...")
        plot23, = ax2.plot(np.multiply(xscale, libKeldysh.FieldToIntensity(Efield1_Keldysh)), PeakWpi_Atomic_Keldysh, '-.k', label=r'$w_{\mathrm{PI}}^{\mathrm{at-anal}}$') #, $m_{\mathrm{eff}}=$'+str(np.round(meff_Keldysh,3))+'.')
        plot2    += [plot23]
    #}}}
    ## ZHUKOV
#    PathPrefix="/home/hilase/Documents/spp-extended-theory/spp_extended_theory/Keldysh/"
    PeakWpi_Zhukov = libKeldyshZhukov.VZ_generateWpiTables(Efield1_Keldysh, 0., pulseWavelength, 800e-9, 0., 0., Egap_SI, meff_Keldysh, Ntotal) #we neglect optical index in atomic theory

    try: #we plot as function of external field
        plot23, = ax2.plot(np.multiply(xscale, libKeldysh.FieldToIntensity(Efield1_Keldysh[np.nonzero(PeakWpi_Zhukov)])), PeakWpi_Zhukov[np.nonzero(PeakWpi_Zhukov)], ':k', label=r'$w_{\mathrm{PI}}^{\mathrm{at-num}}$') #, $m_{\mathrm{eff}}=$'+str(np.round(meff_Waroquiers_LDA,3))+'.') #X axis: external intensity. #Y axis: internal (compressed) intensity.
        plot2    += [plot23]

        if(EnableInsetMagnification):
            plot232, = ax22.plot(np.multiply(xscale, libKeldysh.FieldToIntensity(Efield1_Keldysh[np.nonzero(PeakWpi_Zhukov)])), PeakWpi_Zhukov[np.nonzero(PeakWpi_Zhukov)], ':k', label=r'$w_{\mathrm{PI}}^{\mathrm{at-num}}$')
    except:
        logger.warning(Header+"** Warning: no data available from VP Zhukov for this wavelength.")


    ## GULLEY 
    #plot23, = ax2.loglog(libKeldysh.FieldToIntensity(Efield1_Keldysh), PeakWpi_Gulley, '-.', label=r'$w_{\mathrm{PI}}^{\mathrm{Gu.}}$, $m_{\mathrm{eff}}=$'+str(np.round(meff_Keldysh,3))+'.')
    #plot2    += [plot23]

    ## TUNNELING SECTION

    Egap_SI_Tunnel = Egap_SI #direct band gap for LDA. Indirect requires phonons (or local field effects)

    EnableStarkEffect=False
    if(pulseWavelength  == 800e-9):
        EffectiveMassList=[meff_Waroquiers_LDA]
    elif(pulseWavelength==1600e-9):
        EffectiveMassList=[meff_Waroquiers_LDA]
    elif(pulseWavelength==3200e-9):
        EffectiveMassList=[meff_Waroquiers_LDA]
    else:
        EffectiveMassList=[meff_Waroquiers_LDA]
    #wPI_min={800e-9: 1E36, 1600e-9: 1E35, 2200e-9: 1E30, 2600e-9: 1E30, 3200e-9: 1E34} #desired
    wPI_min={407.58e-9: 1E35, 483e-9: 1E35, 800e-9: 1E36, 1030e-9: 1E36, 1600e-9: 1E35, 2200e-9: 1E30, 2600e-9: 1E30, 3200e-9: 1E34} #tweaked for log
    wPI_max={407.58e-9: 1E45, 483e-9: 1E45, 800e-9: 1E45, 1030e-9: 1E45, 1600e-9: 1E45, 2200e-9: 1E45, 2600e-9: 1E45, 3200e-9: 1E45}


    for meff_Tunnel in EffectiveMassList: #meff_Waroquiers_LDA, meff_VeryHighField]:
        #{{{ #NOTE Stark effect is enabled here, not for the wavelength approach.
        PeakWpi_Tunneling   = np.clip(libKeldysh.KeldyshTunnelingLimit(Egap_SI_Tunnel, meff_Tunnel, pulseWavelength, Efield1_Keldysh*np.sqrt(OpticalIndex), EnableStarkEffect, 1.), 1E30,2E46)
        #NOTE: tunneling barrier does not change as function of refractive index. Hence, we remove the optical index in case of tunneling limit.
        plot23, = ax2.plot(np.multiply(xscale, libKeldysh.FieldToIntensity(Efield1_Keldysh)), PeakWpi_Tunneling[np.nonzero(PeakWpi_Tunneling)], '--k', label=r'$w_{\mathrm{PI}}^{\mathrm{tunnel}}$') #, $m_{\mathrm{eff}}=$'+str(np.round(meff_Tunnel,3))+'.')
        plot2    += [plot23]

        if(EnableInsetMagnification):
            plot232, = ax22.plot(np.multiply(xscale, libKeldysh.FieldToIntensity(Efield1_Keldysh)), PeakWpi_Tunneling[np.nonzero(PeakWpi_Tunneling)], '--k', label=r'$w_{\mathrm{PI}}^{\mathrm{tunnel}}$')

        ## TUNNELING WITH INDIRECT BAND GAP.
        #PeakWpi_Tunneling_i   = np.clip(KeldyshTunnelingLimit(Egap_SI_indirect, meff_Tunnel, pulseWavelength, Efield1_Keldysh, EnableStarkEffect, 1.), 1E30,1E46)
        #plot24, = ax2.loglog(FieldToIntensity(Efield1_Keldysh), PeakWpi_Tunneling_i, label=r'$w_{\mathrm{TI}}^{i}$, $m_{\mathrm{eff}}=$'+str(np.round(meff_Tunnel,3))+'.')
        #plot2    += [plot24]
    #}}}

    #plt.ylim((wPI_min[pulseWavelength],1E45)) #we adjust min of the wPI as function of wavelength
    #plt.xlim(x_min[pulseWavelength])
    ax1.set_xlim(x_min[pulseWavelength])
    ax1.set_ylim([5E20, 1E30])
    ax2.set_xlim(x_min[pulseWavelength])
    ax2.set_ylim([wPI_min[pulseWavelength],wPI_max[pulseWavelength]]) #this removes subticks without asking
    if(EnableInsetMagnification):
        ax22.set_xlim([7.3E10,2E11])
        ax22.set_ylim([3.1E37,2E39])

    # Interpolation of a joint law of photoionization
    #PeakWpi_Effective=np.power(np.power(PeakWpi_Tunneling,-1)+np.power(PeakWpi_Gruzdev,-1),-1)

    #plot23, = ax2.loglog(FieldToIntensity(Efield1_Keldysh), PeakWpi_Effective, 'k-', linewidth=2, label=r'$w_{\mathrm{eff}}$')
    #plot2    += [plot23]

    labels2 = [l.get_label() for l in plot2]
    ax2.legend(plot2, labels2, loc='best', fancybox=True, framealpha=1, fontsize=OneSizeFitsAll, ncol=ncols)
    #ax12.set_yscale('log')
    #ax12.legend(loc='best', fancybox=True, framealpha=0)
    ax2.set_xlabel(r'Intensity (W/cm$^2$)')
    ax2.set_ylabel(r'$w_{\mathrm{PI}}$ (m$^{-3}$ s$^{-1}$)')

    ax1.set_xscale('log')
    ax1.set_yscale('log')

    # Setting the major ticks first (see https://stackoverflow.com/questions/44078409/matplotlib-semi-log-plot-minor-tick-marks-are-gone-when-range-is-large)
    locmaj = mp.ticker.LogLocator(base=10,numticks=100)
    ax1.xaxis.set_major_locator(locmaj)
    locmaj = mp.ticker.LogLocator(base=10,numticks=100)
    ax1.yaxis.set_major_locator(locmaj)

    # Setting the minor ticks then
    # TODO: solution for minor ticks is here: https://stackoverflow.com/questions/44078409/matplotlib-semi-log-plot-minor-tick-marks-are-gone-when-range-is-large
    subdiv=np.arange(0.,10.,1)
    locmin = mp.ticker.LogLocator(base=10.0,subs=subdiv,numticks=1000)
    ax1.xaxis.set_minor_locator(locmin)
    ax1.xaxis.set_minor_formatter(mp.ticker.NullFormatter())
    locmin = mp.ticker.LogLocator(base=10.0,subs=subdiv,numticks=1000)
    ax1.yaxis.set_minor_locator(locmin)
    ax1.yaxis.set_minor_formatter(mp.ticker.NullFormatter())

    #ax1.xaxis.set_minor_locator(mp.ticker.LogLocator(base=10.0, subs=None, numticks=10, numdecs=None))
    #ax1.yaxis.set_minor_locator(mp.ticker.LogLocator(base=10.0, subs=None, numticks=10, numdecs=None))
    #ax1.minorticks_off()

    ax2.set_xscale('log') #compared to using loglog(), the xtics and ytics are better for publication.
    ax2.set_yscale('log')
    if(EnableInsetMagnification):
        ax22.set_xscale('log')
        ax22.set_yscale('log')

    # Setting major ticks locators
    locmaj = mp.ticker.LogLocator(base=10,numticks=100)
    ax2.xaxis.set_major_locator(locmaj)
    locmaj = mp.ticker.LogLocator(base=10,numticks=100)
    ax2.yaxis.set_major_locator(locmaj)

    if(not SidePanels): #we remove the ticks and marks from top figure
        ax1.xaxis.set_major_formatter(mp.ticker.NullFormatter())
        ax1.xaxis.set_minor_formatter(mp.ticker.NullFormatter())
        ax1.set_xlabel("")
    if(EnableInsetMagnification):
        locmaj = mp.ticker.LogLocator(base=10,numticks=2)
        ax22.xaxis.set_major_locator(locmaj)
        locmaj = mp.ticker.LogLocator(base=10,numticks=2)
        ax22.yaxis.set_major_locator(locmaj)

    # Setting minor ticks locators
    locmin = mp.ticker.LogLocator(base=10.0,subs=subdiv,numticks=1000)
    ax2.xaxis.set_minor_locator(locmin)
    ax2.xaxis.set_minor_formatter(mp.ticker.NullFormatter())
    locmin = mp.ticker.LogLocator(base=10.0,subs=subdiv,numticks=1000)
    ax2.yaxis.set_minor_locator(locmin)
    ax2.yaxis.set_minor_formatter(mp.ticker.NullFormatter())

    if(EnableInsetMagnification):
        locmin = mp.ticker.LogLocator(base=10.0,subs=subdiv,numticks=2)
        ax22.xaxis.set_minor_locator(locmin)
        #ax22.xaxis.set_minor_formatter(mp.ticker.NullFormatter())
        locmin = mp.ticker.LogLocator(base=10.0,subs=subdiv,numticks=2)
        ax22.yaxis.set_minor_locator(locmin)
        #ax22.yaxis.set_minor_formatter(mp.ticker.NullFormatter())

    #ax2.xaxis.set_minor_locator(mp.ticker.LogLocator(base=10.0, subs=None, numticks=10, numdecs=None))
    #ax2.yaxis.set_minor_locator(mp.ticker.LogLocator(base=10.0, subs=None, numticks=10, numdecs=None))
    #ax2.minorticks_off()

    #print(ax2.yaxis.minorTicks.count)
    #ax2.yaxis.set_major_formatter()

    #ax2.yaxis.set_major_locator(mp.ticker.LogLocator(numticks=10))
    #ax2.legend(plot1, labels1, loc='best', fancybox=True, framealpha=0)
    #ax2.axes.get_xaxis().set_visible(False)
    #ax2.axes.get_yaxis().set_visible(False)
    #else:
    #plt.tight_layout()

    ## This should be the inset that we would put in Fig. 3a.
    # Previously, it was prepared using GNUplot file plotFieldWpiPeaks.gnu.

    #if(PlotVirtualExcitation): #{{{ #NOTE: no longer used, since it has no meaning. Gauge dependent quantity.

    #ax3 = plt.subplot(plotC)
    #plot1 = [] #np.empty((0,0))
    #for pulse_value in pulses_uniq: #{{{  
    ##plt.title('Virtual exc. factor')
    ## select the corresponding wavelength in results
    #PlottingSinglePulses_pulses = libDatabase.FilterDatabase(PlottingSinglePulses, pulse_value, IndexTable['pulses1'])
    #PlottingSinglePulses_pulses = libDatabase.FilterDatabase(PlottingSinglePulses_pulses, pulse_value, IndexTable['pulses2']) #keeps only pulses of the same duration.

    #runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, ElectronEnergy_max, Energy_max, W_pi, PulseType, RiseTime, FinalTime = libPlotting.ExtractData(PlottingSinglePulses_pulses)

    #print Header+"** Info: saving the data in CSV file for fitting. "
    #filename = "EffectOfPulseIntensity-Duration-"+str(pulseWavelength*1E9)+"nm-tau"+str(int(pulse_value*1E15))+"fs"
    #print Header+"** Info: saving in "+filename+".csv"
    ##libDatabase.ExportToTxt(PlottingSinglePulses_pulses, filename+".csv")
    #lists = sorted(itertools.izip(*[libKeldysh.FieldToIntensity(Efield1), irreversibility, pulses1]))
    #Intensity_sorted, irreversibility_sorted, pulses1_sorted = list(itertools.izip(*lists))

    #plot_field, = ax2.plot(1E-4*Intensity_sorted, irreversibility_sorted, Colors[pulse_value]+'-o', label=r'$\zeta$, $\tau=$'+str(round(pulse_value*1E15))+' fs')
    #plot1 += [plot_field]

    ##}}}
    #plt.xlabel(r'Peak intensity (W/cm$^2$)')
    #plt.ylabel(r'Virtual exc. $\zeta$')
    #plt.ylim((0,1))
    #ax3.set_xscale('log')
    #ax3.grid()

    #ax32 = ax3.twinx()
    #for pulse_value in pulses_uniq:
    ## select the corresponding wavelength in results
    #PlottingSinglePulses_pulses = libDatabase.FilterDatabase(PlottingSinglePulses, pulse_value, IndexTable['pulses1']) #keeps only 800 nm wavelenghts.

    #runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, ElectronEnergy_max, Energy_max, W_pi, PulseType, RiseTime, FinalTime = libPlotting.ExtractData(PlottingSinglePulses_pulses)

    #plot_field, = ax32.plot(FieldToIntensity(Efield1), gammaKeldysh(Egap_SI, meff_Waroquiers_LDA, Efield1, wavelengths1), 'k-', label=r'$\gamma$')
    ##plot1 += [plot_field]

    #plt.ylabel(r'Adiabadicity $\gamma$')
    ##plt.ylim((0,1))
    #ax32.set_yscale('log')

    ## ax2.set_yscale("log")

    #labels3 = [l.get_label() for l in plot3]
    ##ax22.legend(plot1, labels1, loc='center right', fancybox=True, framealpha=0)
    ##ax1.legend(loc=1, borderaxespad=0.)
    ## ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ##ax22.legend(loc=2, borderaxespad=0.)
    #plt.tight_layout()
    ##}}}


    def FitDrudeModel(n_exc, nu):
        #omega       = 2.*np.pi*c/wavelength
        omega_p2    = np.multiply(n_exc, e**2 / m_e / meff / epsilon_0)
        Epsilon = 1.  - np.divide(omega_p2, np.power(omega, 2)) * np.divide(1.,1. + np.divide(1.j * nu, omega)) #conduction band absorption, no interband contribution
        return Epsilon

    # d e f FitDrudeModelReal(n_exc, nu):
    #Epsilon = FitDrudeModel(n_exc, nu)
    ##1.  - np.divide(omega_p2, np.power(omega, 2)) * np.divide(1.,1. + np.divide(1.j * nu, omega)) #conduction band absorption, no interband contribution
    ## the coupling/interference between interband and intraband motions is accounted for in the complete dielectric function!
    #DrudeEnergy          = SuperCellVolume_SI / e * epsilon_0 * omega * Epsilon.imag * np.multiply(Efield_X_SI_interp_Ne, Efield_X_SI_interp_Ne)
    #TotalDrudeEnergy     = np.multiply(np.cumsum(DrudeEnergy), dt_deriv[0])+ExcessEnergy

    #return TotalDrudeEnergy


    if(PlotEnergy):
        ax3 = plt.subplot(plotC)
        #plt.title('Virtual exc. factor')

        plot3 = []
        for pulse_value in pulses_uniq: #{{{
            # select the corresponding wavelength in results
            PlottingSinglePulses_pulses = libDatabase.FilterDatabase(PlottingSinglePulses, pulse_value, libPlotting.IndexTable['pulses1'])
            PlottingSinglePulses_pulses = libDatabase.FilterDatabase(PlottingSinglePulses_pulses, pulse_value, libPlotting.IndexTable['pulses2']) #keeps only pulses of the same duration.

            runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingTable)

            if(UseNexcAfterPulse):
                del Ne_final
                Ne_final = Ne_afterpulse #NOTE: dirty and valid only until next call of ExtractData

            lists = sorted(zip(*[libKeldysh.FieldToIntensity(Efield1), ElectronEnergy_max, pulses1, Ne_final]))
            Intensity_sorted, Energy_final_sorted, pulses1_sorted, Ne_final_sorted = list(zip(*lists))
            #print Header+"Ready to plot: "

            plot_field, = ax3.plot(np.multiply(xscale, Intensity_sorted), Energy_final_sorted, Colors[pulse_value]+'-'+Symbols[pulse_value], label=r'$\tau=$'+str(int(pulse_value*1E15))+" fs")
            plot3 += [plot_field]
            ## TODO: can we adjust collision frequency on Drude here?

            #try:
            #FitDrudeModel_parameters = curve_fit(FitDrudeModelReal, Ne_final_sorted, Energy_final_sorted, bounds=(1E13,1E16))
            # We have to minimize ElectronEnergy -
            #print Header+"FITTING DRUDE MODEL."
            #print np.divide(1.,FitDrudeModel_parameters)
            #nu_t = FitDrudeModel_parameters[0]

            #TODO: still have to use epsilon_Drude to compute Energy.
            #plot_field, = ax3.plot(Intensity_sorted, FitDrudeModel(Ne_final_sorted,nu_t), Colors[pulse_value]+'-.', label=r'Drude, $\tau=$'+str(int(pulse_value*1E15))+r" fs, $\nu^{-1}=$"+str(round(nu_t,2)))
            #plot3 += [plot_field]
            #except:
            #print Header+"** WARNING: auto-fit of collision frequency failed."
            #nu_t = []


        #}}}
        labels3 = [l.get_label() for l in plot3]
        ax3.legend(plot3, labels3, loc='best', fancybox=True, framealpha=1, fontsize=OneSizeFitsAll)
        plt.xlabel(r'Peak intensity (W/cm$^2$)')
        plt.ylabel(r'Energy $\xi_{\mathrm{el}}^{\mathrm{TDDFT}}$ (eV)')
        #plt.xlim((0,1))
        #locmin = mp.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1, 2, 3,4,5,6,7,8,9,10))
        #locmin = mp.ticker.LogLocator(base=10.0,numticks=4)
        #ax1.get_xaxis().set_major_locator(mp.ticker.LogLocator(base=10.0,numticks=4))
        #ax1.get_xaxis().set_minor_locator(mp.ticker.LogLocator(base=10.0,numticks=50))
        #ax1.xaxis.set_minor_locator(mp.ticker.AutoMinorLocator())
        #ax1.get_xaxis().set_minor_formatter(mp.ticker.NullFormatter())
        #ax1.get_xaxis().set_minor_formatter(mp.ticker.ScalarFormatter())


        #ax2.xaxis.set_minor_locator(locmin)
        #ax2.xaxis.set_minor_formatter(mp.ticker.NullFormatter())

        #ax3.xaxis.set_minor_locator(locmin)
        #ax3.xaxis.set_minor_formatter(mp.ticker.NullFormatter())
        #ax3.yaxis.set_major_locator(locmin)
        #ax3.yaxis.set_minor_formatter(mp.ticker.NullFormatter())

        ax3.set_xscale('log')
        #ax3.set_yscale('symlog', linthreshy=1e-5)
        ax3.set_yscale('log')

        plt.tight_layout()
        #ax3.grid()
        #ax32 = ax3.twinx()
        #for pulse_value in pulses_uniq:
        # select the corresponding wavelength in results
        #PlottingSinglePulses_pulses = libDatabase.FilterDatabase(PlottingSinglePulses, pulse_value, IndexTable['pulses1']) #keeps only 800 nm wavelenghts.

        #runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, ElectronEnergy_max, Energy_max, W_pi, PulseType, RiseTime = libPlotting.ExtractData(PlottingSinglePulses_pulses)

        #plt.plot(FieldToIntensity(Efield1), gammaKeldysh(2.56*e, 0.228, Efield1, wavelengths1), '--', label=r'Adiabaticity, $\tau=$'+str(round(pulse_value*1E15))+' fs')
        #plt.ylabel('Adiabadicity factor')
        #plt.ylim((0,1))
        #ax32.set_yscale('log')

        # ax2.set_yscale("log")
        #ax3.legend(loc=1, borderaxespad=0.)
        # ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        #ax32.legend(loc=2, borderaxespad=0.)

    # Saving the plot
    plt.tight_layout()

    if(PlotSaturation):
        Filename_Sat = "-WithSaturation"
    else:
        Filename_Sat = ""

    ExportName = "EffectOfPulseIntensity-Duration-"+str(pulseWavelength*1E9)+"nm"+MoveIonFile+FileName_opt+Filename_Sat
    plt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.eps')
    plt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.png')
    logger.info(Header+"** Saved: 'Metaplot-SinglePulse-SingleColor-'"+ExportName+".png.")
    if(ShowFinalPlot==1):
        plt.show()


    del PlottingSinglePulses, wavelength_uniq, CEP1_uniq
    del ax2
    if(PlotEnergy):
        del ax3
#}}}

def plotEffectOfIntensity_PulseDuration_Metals(PlottingTable, pulseWavelength=800e-9, CEP0=0., PulseType=1, RiseCycles=2, PlotWpiAxis=False, ShowFinalPlot=0, area=1., OpticalIndex=1.0, RatioExternal=1., MoveIon=0): #{{{
  PlotVirtualExcitation = False #Not advised! This is not a real observable. 
  PlotEnergy = True #PRL: False. #PRB: True.
  ShowPulseDurations=True
  PlotFrozenHxc=False
  UseFinalInstantData=True #True: we plot results at final time. False: we plot results right after the pulse. 
  ShowPoynting        = True
  EnableInducedFields = False #True: useful for plasmonics and spectroscopy. False: useful for energy absorption. 
  MinPulseDuration=9E-15
  meff = 1. 
  
  Header="[plotFinalQuantities: plotEffectOfIntensity_PulseDuration_Metals()] "
  
  if(MoveIon): 
      MoveIonTitle=", w/ phonons."
      MoveIonFile="MoveIon"
  else: 
      MoveIonTitle=""
      MoveIonFile=""
  
  if(PlotFrozenHxc):
      FreezeHXC_name="-FrozenHXC-"
  else: 
      FreezeHXC_name="" 
  
  
  print("====== Plot: SINGLE PULSE, SAME WAVELENGTHS, SAME POLARIZATIONS: final density as function of field and pulse duration ======")
  #Filter: double pulses, same wavelengths, same polarization
  PlottingSinglePulses = PlottingTable #get the full database
  PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, libPlotting.IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
  PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, libPlotting.IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0., libPlotting.IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, CEP0, libPlotting.IndexTable['CEP1']) #removes fields where CEP != 0. #NOTE: Can be done for Si-LDA.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 1., libPlotting.IndexTable['TheoryLevel']) #removes runs where TheoryLevel != dft.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, PulseType, libPlotting.IndexTable['PulseType']) #keeps only pulses of type SinSquared (1)
  #PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, RiseCycles, libPlotting.IndexTable['RiseTime']) #keeps only pulses of a selected RiseCycles
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, pulseWavelength, libPlotting.IndexTable['wavelengths1'])
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, MoveIon, libPlotting.IndexTable['MoveIon']) #filters the enabled or disabled movement of ions.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0., libPlotting.IndexTable['GaugeFieldAmplitude']) #filters the probe pulses
  
  PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, MinPulseDuration, libPlotting.IndexTable['pulses1'])
  
  #wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = libPlotting.CountTheNumberOfLaserParameters(PlottingSinglePulses)
  logger.info((Header+"** Info: Number of entries: "+str(len(PlottingSinglePulses))))

  #print("Info: available pulse durations:\n"+str(pulses_uniq))
  
  Symbols={20e-15: '+', 30e-15: 'o', 60e-15: '^'}
  Colors= {20e-15: 'g', 30e-15: 'r', 60e-15: 'b'}
  
  plt.figure(figsize=(10,6))
  ax1 = plt.subplot(211)
  ax2 = plt.subplot(212)
  ax22= ax2.twinx()
  plot1 = []; plot2=[]
  
  WithInducedFields_txt=""
  if(EnableInducedFields):
      PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 1, libPlotting.IndexTable['FieldInduced'])
      WithInducedFields_txt="-WithInducedFields"
  else:
      PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0, libPlotting.IndexTable['FieldInduced'])
      WithInducedFields_txt="-NoInducedFields"
  
  # Loading available experimental data
  
  if(pulseWavelength == 800e-9): 
      # open the experimental file
      ElectronEnergy_experimental_sorted = np.loadtxt("Shugaev-2018-Exp-FluenceReflectivity.csv", delimiter="\t", skiprows=1)
      logger.info((ElectronEnergy_experimental_sorted[0,:]))
  
  
  logger.info((Header+"** Info: plotting full TDDFT..."))
  if(not PlotFrozenHxc):
    PlottingSinglePulses_Normal = libDatabase.FilterDatabase(PlottingSinglePulses, 0., libPlotting.IndexTable['FreezeHXC']) #removes fields where FreezeHXC != 0.
    

    wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = libPlotting.CountTheNumberOfLaserParameters(PlottingSinglePulses)
    logger.info((Header+"** Info: Number of entries: "+str(len(PlottingSinglePulses))))
    logger.info(("Info: available pulse durations:\n"+str(pulses_uniq)))

    
    logger.info((Header+"** Info: Number of entries after FreezeHXC: "+str(len(PlottingSinglePulses_Normal))))
    for pulse_value in pulses_uniq: #{{{
            # select the corresponding wavelength in results
            PlottingSinglePulses_pulses = libDatabase.FilterDatabase(PlottingSinglePulses_Normal, pulse_value, libPlotting.IndexTable['pulses1'])
            logger.info((Header+"** Info: Number of entries: "+str(len(PlottingSinglePulses_pulses))))

            #PlottingSinglePulses_pulses = libDatabase.FilterDatabase(PlottingSinglePulses_pulses, pulse_value, IndexTable['pulses2']) #keeps only pulses of the same duration.
            
            wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = libPlotting.CountTheNumberOfLaserParameters(PlottingSinglePulses_pulses)
            
            runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingTable)

           
            ## build the internal laser peak fluence here
            intensity1 = 0.5*c*epsilon_0*Efield1**2 #W/m2
            fluence1   = intensity1*pulses1 #J/m2 #NOTE: exact while sin² envelope is used. 
            
            #RatioExternal=24 #given by Yabana et al., PRB 2012. # TODO: detail!!!
            if(EnableInducedFields):
                intensity1_ext = 0.5*c*epsilon_0*(Efield1/RatioExternal)**2 #W/m2 # WARNING : may be false!
            #intensity1_ext = 0.5*c*epsilon_0*(RatioExternal*Efield1)**2 #W/m2 
            else:
                intensity1_ext = 0.5*c*epsilon_0*(Efield1)**2 #W/m2 
            fluence_ext1 = intensity1_ext*pulses1 #J/m2 #still exact
            
            ## build the external laser peak fluence here
            
            ## build the laser peak fluence here
            incident_energy1 = fluence1 * area / e #WARNING: to be verified. Should be in eV.
            incident_energy1_ext = fluence_ext1 * area / e
            
            if(UseFinalInstantData):
                ElectronEnergy_use    = ElectronEnergy_max
                EigenvaluesEnergy_use = EigenvaluesEnergy_max
                CaptionText           = r'$(t_{\mathrm{end}})$'
            else: 
                ElectronEnergy_use    = ElectronEnergy_afterpulse
                EigenvaluesEnergy_use = EigenvaluesEnergy_afterpulse
                CaptionText           = r'$(t \geq \tau)$'
            
            ## absorbed energy relative to laser pulse energy
            relative_absorbed_energy1 = np.divide(ElectronEnergy_use, incident_energy1_ext) #eV over eV
            #relative_absorbed_energy1 = np.divide(ElectronEnergy_max, incident_energy1_ext) #eV over eV
            relative_absorbed_energy_eigenvalues1 = np.divide(EigenvaluesEnergy_use, incident_energy1_ext) #eV over eV
            
            #lists = sorted(itertools.izip(*[libKeldysh.FieldToIntensity(Efield1), fluence1, fluence_ext1, ElectronEnergy_max, pulses1, Ne_final]))        
            lists = sorted(zip(*[pulses1, libKeldysh.FieldToIntensity(Efield1), fluence_ext1, ElectronEnergy_use,  EigenvaluesEnergy_use, EnergyPoynting_SI_sum_final, relative_absorbed_energy1, relative_absorbed_energy_eigenvalues1]))
            
            ## WARNING : Does this create negative values? 
            pulses1_sorted, Intensity_sorted, fluence_ext1_sorted, ElectronEnergy_use_sorted, EigenvaluesEnergy_use_sorted, EnergyPoynting_SI_sum_final_sorted, relative_absorbed_energy1_sorted, relative_absorbed_energy_eigenvalues1_sorted = list(zip(*lists))
            #print Header+"Ready to plot: "
            
            plot_field, = ax1.plot(np.multiply(xscale, Intensity_sorted), ElectronEnergy_use_sorted, Colors[pulse_value]+'^', label=r'abs.'+CaptionText+r', $\tau=$'+str(int(pulse_value*1E15))+" fs")
            plot1 += [plot_field]
            
            logger.debug((Header+"** Debug: contents of Energy Poynting table..."))
            logger.debug(EnergyPoynting_SI_sum_final_sorted)
            
            if(ShowPoynting):
                plot_field, = ax1.plot(np.multiply(xscale, Intensity_sorted), EnergyPoynting_SI_sum_final_sorted, Colors[pulse_value]+'-', label=r"$\int J(t').E(t') dt' $"+r', $\tau=$'+str(int(pulse_value*1E15))+" fs")
                plot1 += [plot_field]
            
            #plot_field, = ax1.plot(np.multiply(xscale, Intensity_sorted), ElectronEnergy_afterpulse_sorted, Colors[pulse_value]+'-o', label=r'abs.$(t_{\mathrm{end}})$, $\tau=$'+str(int(pulse_value*1E15))+" fs")
            #plot1 += [plot_field]
            
            plot_field, = ax1.plot(np.multiply(xscale, Intensity_sorted), EigenvaluesEnergy_use_sorted, Colors[pulse_value]+'--^', label=r'eig.'+CaptionText+r', $\tau=$'+str(int(pulse_value*1E15))+" fs")
            plot1 += [plot_field]
            
            
            #plot_field, = ax1.plot(np.multiply(xscale, Intensity_sorted), EigenvaluesEnergy_afterpulse_sorted, Colors[pulse_value]+'-.', label=r'eig. $t>\tau$, $\tau=$'+str(int(pulse_value*1E15))+" fs")
            #plot1 += [plot_field]
            
            #plot_field, = ax2.plot(np.multiply(xscale, fluence1_sorted), relative_absorbed_energy1_sorted, Colors[pulse_value]+'--'+Symbols[pulse_value], label=r'$\xi_{el}(F_{int})$, $\tau=$'+str(int(pulse_value*1E15))+" fs")
            #plot2 += [plot_field]
            
            plot_field, = ax2.plot(np.multiply(xscale, fluence_ext1_sorted), relative_absorbed_energy1_sorted, Colors[pulse_value]+'-'+Symbols[pulse_value], label=r'$\xi_{el}$'+CaptionText+r'$[F_{ext}]$, $\tau=$'+str(int(pulse_value*1E15))+" fs")
            plot2 += [plot_field]
            
            if(pulseWavelength == 800e-9 and pulse_value==60e-15): 
                plot_field, = ax22.plot(ElectronEnergy_experimental_sorted[:,0], 1.-ElectronEnergy_experimental_sorted[:,1], Symbols[pulse_value], label=r'exp. abs.'+CaptionText+r', $\tau=$'+str(int(pulse_value*1E15))+" fs")
                plot2 += [plot_field]
            
            #plot_field, = ax2.plot(np.multiply(xscale, fluence_ext1_sorted), relative_absorbed_energy_eigenvalues1_sorted, Colors[pulse_value]+'--'+Symbols[pulse_value], label=r'eig., $\xi_{eig}(F_{ext})$, $\tau=$'+str(int(pulse_value*1E15))+" fs")
            #plot2 += [plot_field]
            

            del lists
  #}}}
  #if(PlotFrozenHxc):
  else:
    ## We can now filter and plot the FreezeHXC results

    logger.info((Header, "** Info: plotting effect of Hxc..."))
    PlottingSinglePulses_HXC = libDatabase.FilterDatabase(PlottingSinglePulses, 1., libPlotting.IndexTable['FreezeHXC']) #keep only fields where FreezeHXC != 0.
    
    wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = libPlotting.CountTheNumberOfLaserParameters(PlottingSinglePulses)
    logger.info((Header+"** Info: Number of entries: "+str(len(PlottingSinglePulses))))
    logger.info(("Info: available pulse durations:\n"+str(pulses_uniq)))
    
    for pulse_value in pulses_uniq: #{{{
            # select the corresponding wavelength in results
            PlottingSinglePulses_pulses = libDatabase.FilterDatabase(PlottingSinglePulses_HXC, pulse_value, libPlotting.IndexTable['pulses1'])
            #PlottingSinglePulses_pulses = libDatabase.FilterDatabase(PlottingSinglePulses_pulses, pulse_value, libPlotting.IndexTable['pulses2']) #keeps only pulses of the same duration.

            wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = libPlotting.CountTheNumberOfLaserParameters(PlottingSinglePulses_pulses)
            
            if(len(Efield1_uniq)==0): 
                logger.warning((Header, "** Warning: no data available for this run."))
                continue

            runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingSinglePulses_pulses)
            
            ## build the internal laser peak fluence here
            intensity1 = 0.5*c*epsilon_0*Efield1**2 #W/m2
            fluence1   = intensity1*pulses1 #J/m2 #NOTE: exact while sin² envelope is used. 
            
            #RatioExternal=24 #given by Yabana et al., PRB 2012. # TODO: detail!!!
            if(EnableInducedFields):
                intensity1_ext = 0.5*c*epsilon_0*(Efield1/RatioExternal)**2 #W/m2 #ERROR?
            #intensity1_ext = 0.5*c*epsilon_0*(RatioExternal*Efield1)**2 #W/m2
            else:
                intensity1_ext = 0.5*c*epsilon_0*(Efield1)**2
            fluence_ext1 = intensity1_ext*pulses1 #J/m2 #still exact
            
            ## build the external laser peak fluence here
            
            ## build the laser peak fluence here
            incident_energy1 = fluence1 * area / e #WARNING: to be verified. Should be in eV.
            incident_energy1_ext = fluence_ext1 * area / e
            
            if(UseFinalInstantData):
                ElectronEnergy_use    = ElectronEnergy_max
                EigenvaluesEnergy_use = EigenvaluesEnergy_max
                CaptionText           = r'$(t_{\mathrm{end}})$'
            else: 
                ElectronEnergy_use    = ElectronEnergy_afterpulse
                EigenvaluesEnergy_use = EigenvaluesEnergy_afterpulse
                CaptionText           = r'$(t > \tau)$'
        
            ## absorbed energy relative to laser pulse energy
            relative_absorbed_energy1 = np.divide(ElectronEnergy_use, incident_energy1_ext) #eV over eV
            #relative_absorbed_energy1 = np.divide(ElectronEnergy_max, incident_energy1_ext) #eV over eV
            relative_absorbed_energy_eigenvalues1 = np.divide(EigenvaluesEnergy_use, incident_energy1_ext) #eV over eV
            
            #lists = sorted(itertools.izip(*[libKeldysh.FieldToIntensity(Efield1), fluence1, fluence_ext1, ElectronEnergy_max, pulses1, Ne_final]))        
            lists = sorted(zip(*[pulses1, libKeldysh.FieldToIntensity(Efield1), fluence_ext1, ElectronEnergy_use,  EigenvaluesEnergy_use, relative_absorbed_energy1, relative_absorbed_energy_eigenvalues1]))
            
            ## WARNING : Does this create negative values? 
            pulses1_sorted, Intensity_sorted, fluence_ext1_sorted, ElectronEnergy_use_sorted, EigenvaluesEnergy_use_sorted, relative_absorbed_energy1_sorted, relative_absorbed_energy_eigenvalues1_sorted = list(zip(*lists))
            #print Header+"Ready to plot: "
            
            plot_field, = ax1.plot(np.multiply(xscale, Intensity_sorted), ElectronEnergy_use_sorted, Colors[pulse_value]+'-', label=r'abs.'+CaptionText+r', $H_{xc}$ frozen, $\tau=$'+str(int(pulse_value*1E15))+r" fs")
            plot1 += [plot_field]
            
            #plot_field, = ax1.plot(np.multiply(xscale, Intensity_sorted), ElectronEnergy_afterpulse_sorted, Colors[pulse_value]+'-o', label=r'abs.$(t_{\mathrm{end}})$, $\tau=$'+str(int(pulse_value*1E15))+" fs")
            #plot1 += [plot_field]
            
            plot_field, = ax1.plot(np.multiply(xscale, Intensity_sorted), EigenvaluesEnergy_use_sorted, Colors[pulse_value]+'^', label=r'eig.'+CaptionText+r', $H_{xc}$ frozen, $\tau=$'+str(int(pulse_value*1E15))+r" fs")
            plot1 += [plot_field]
            
            #plot_field, = ax1.plot(np.multiply(xscale, Intensity_sorted), EigenvaluesEnergy_afterpulse_sorted, Colors[pulse_value]+'-.', label=r'eig. $t>\tau$, $\tau=$'+str(int(pulse_value*1E15))+" fs")
            #plot1 += [plot_field]
            
            #plot_field, = ax2.plot(np.multiply(xscale, fluence1_sorted), relative_absorbed_energy1_sorted, Colors[pulse_value]+'--'+Symbols[pulse_value], label=r'$\xi_{el}(F_{int})$, $\tau=$'+str(int(pulse_value*1E15))+" fs")
            #plot2 += [plot_field]
            
            plot_field, = ax2.plot(np.multiply(xscale, fluence_ext1_sorted), relative_absorbed_energy1_sorted, Colors[pulse_value]+'-'+Symbols[pulse_value], label=r'$\xi_{el}$'+CaptionText+r'$[F_{ext}]$, $H_{xc}$ frozen, $\tau=$'+str(int(pulse_value*1E15))+r" fs")
            plot2 += [plot_field]
            
            #plot_field, = ax2.plot(np.multiply(xscale, fluence_ext1_sorted), relative_absorbed_energy_eigenvalues1_sorted, Colors[pulse_value]+'--'+Symbols[pulse_value], label=r'eig., $\xi_{eig}(F_{ext})$, $\tau=$'+str(int(pulse_value*1E15))+" fs")
            #plot2 += [plot_field]
            
            del lists
  #}}}
  
  labels1 = [l.get_label() for l in plot1]
  labels2 = [l.get_label() for l in plot2]
  ax1.legend(plot1, labels1, loc='best', fancybox=True, framealpha=1, fontsize=OneSizeFitsAll)
  ax2.legend(plot2, labels2, loc='best', fancybox=True, framealpha=1, fontsize=OneSizeFitsAll)
  
  ax1.set_xscale('log')
  #ax1.set_yscale('log')
  ax1.grid()
  ax1.set_yscale('symlog', linthreshy=1e-3)
  
  ax1.set_xlabel(r'Peak intensity (W/cm$^2$)')
  ax1.set_ylabel(r'Energy $\xi_{\mathrm{el}}^{\mathrm{TDDFT}}$ (eV)')
  
  ax2.set_xscale('log')
  #ax2.set_xlim([7E8, 2E15])
  #ax2.set_yscale('symlog', linthreshy=1e-10)
  ax2.grid()
  ax2.set_xlabel(r'Peak fluence (J/cm$^2$)')
  ax2.set_ylabel(r'Rel. abs. $\xi_{\mathrm{el}}^{\mathrm{TDDFT}} / \xi_{\mathrm{laser}}$')
  ax22.set_ylabel(r'Exp. $1-R$')
  
  #ax1.yaxis.set_major_locator(mp.ticker.LogLocator(base=10.0, subs='auto', numticks=None))
  ax1.yaxis.set_major_locator(mp.ticker.FixedLocator([-1E4, -1E2, -10, -1, -1E-2, -1E-3, 0, 1E-3, 1E-2, 1E0, 1E2, 1E4])) #base=10.0, subs='auto', numticks=None))
  
  
  plt.tight_layout()
  #print(Header, "Debug: PulseType: ", PulseType)
  ExportName = "EffectOfPulseIntensity-Duration-PulseType"+str(libPlotting.PulseTypeID_rev[PulseType[0]])+"-Wavelength"+str(pulseWavelength*1E9)+"nm"+FreezeHXC_name+MoveIonFile+WithInducedFields_txt
  plt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.eps')
  plt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.png')
  logger.info(Header+"** Saved: 'Metaplot-SinglePulse-SingleColor-'"+ExportName+".png.")
  
  if(ShowFinalPlot==1):
    plt.show()
#}}}
    
## Single-pulse: Final electron density as function of intensity. Each curve is for a different CEP. 
# TODO: current and energy, and harmonic spectrum as function of CEP.
def plotEffectOfIntensity_CEP(PlottingTable, pulseWavelength=800e-9, pulseDuration=20e-15, FilterCEP=False, CEP0=0., PulseShape=3, pulseRiseTime=0, ShowFinalPlot=0, Egap_SI=2.56*e): #{{{
  print("#====== Plot: SINGLE PULSE, EFFECT OF THE CEP, SAME WAVELENGTHS, SAME POLARIZATIONS: final density as function of field ======")
  #Filter: double pulses, same wavelengths, same polarization

  # Simple pulse, single color
  PlottingSinglePulses = PlottingTable #get the full database
  PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, libPlotting.IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
  PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, libPlotting.IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0., libPlotting.IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, pulseWavelength, libPlotting.IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0, libPlotting.IndexTable['FreezeHXC']) #keeps only non-frozen XC.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, PulseShape, libPlotting.IndexTable['PulseType']) #keeps only non-frozen XC.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, pulseDuration, libPlotting.IndexTable['pulses1']) #removes fields where tau != pulseDuration

  if(FilterCEP==True):
      PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, CEP0, libPlotting.IndexTable['CEP1']) #removes fields where tau != pulseDuration
  if(pulseRiseTime != 0): 
    PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, pulseRiseTime, libPlotting.IndexTable['RiseTime']) #removes fields where tau != pulseDuration

  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = libPlotting.CountTheNumberOfLaserParameters(PlottingSinglePulses)
  
  plt.figure(figsize=(SizeX,SizeY))
  ax1 = plt.subplot(111) #111 or plotA if using plotB
  plt.title(r'$n_{\mathrm{exc}}(I)$, $\tau=$'+str(int(pulseDuration*1E15))+' fs.')
  # Number of sweeped parameters: 3 (Efield, Delay, CrossPolarization)
  #plt.scatter(libKeldysh.FieldToIntensity(Efield1), Ne_final, c=Delay, s=(30.+CrossPolarization*50.), label=r'$N_{exc}$[Field]')
  plot1 = [] #np.empty((0,0))
  for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
    PlottingSinglePulses_CEP = libDatabase.FilterDatabase(PlottingSinglePulses, CEP1_value, libPlotting.IndexTable['CEP1']) #select the local CEP
    
    ## prints the current state of selection and stores the values of sweeped parameters

    runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingTable)

    
    #plot_field, = ax1.plot(libKeldysh.FieldToIntensity(Efield1), Ne_final, '-', label=r'$N_{exc}^{final}$, $\phi_{CE}$='+str(round(CEP1_value*180./pi)))
    #lists = sorted(itertools.izip(*[libKeldysh.FieldToIntensity(Efield1), Ne_final]))
    lists = sorted(zip(*[libKeldysh.FieldToIntensity(Efield1), Ne_final, pulses1]))
    Intensity_sorted, Ne_final_sorted, pulses1_sorted = list(zip(*lists))
    #print Header+"Ready to plot: "
    #print "** Info: Sorting data..."
    #print Efield1_sorted, Ne_final_sorted
    plot_field, = ax1.plot(Intensity_sorted, Ne_final_sorted, 'o', label=r'$N_{exc}^{final}$, $\phi_{CE}$='+str(round(CEP1_value*180./pi)))
    logger.info(Header+"plotEffectOfIntensity_CEP(): adding CEP "+str(round(CEP1_value*180./pi))+" to the the plot.")
    logger.info(plot_field)
    logger.debug(Header+"{Fore.RED}** Debug info: "+str(np.shape(plot_field)))
    plot1 += [plot_field] # does cumulate each CEP now. Type: Line2D, not np.array. 
    
    #ax1.legend(loc=4, borderaxespad=0.) #bbox_to_anchor=(1.05, 1), 
  labels1 = [l.get_label() for l in plot1]
  ax1.legend(plot1, labels1, loc=2)
  plt.xlabel(r'Peak intensity (W/cm$^2$)')
  plt.ylabel(r'Final density (m$^{-3}$)')
  #  plt.ylim((1E16,4.*5e28))
  ax1.set_xscale('log')
  ax1.set_yscale('log')
  #plot1_list = plot1.tolist()
    #print Efield1, CEP1_value
    
  #Adding axis of gammaKeldysh, but for all CEPs !!
  #runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, ElectronEnergy_max, Energy_max, W_pi, PulseType, RiseTime = libPlotting.ExtractData(PlottingSinglePulses)
  
  ax12 = ax1.twinx()
  for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
    PlottingSinglePulses_CEP = libDatabase.FilterDatabase(PlottingSinglePulses, CEP1_value, libPlotting.IndexTable['CEP1']) #select the local CEP
    
    ## prints the current state of selection and stores the values of sweeped parameters
    runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingTable)
    
    #plot_field, = ax1.plot(libKeldysh.FieldToIntensity(Efield1), Ne_final, '-', label=r'$N_{exc}^{final}$, $\phi_{CE}$='+str(round(CEP1_value*180./pi)))
    #lists = sorted(itertools.izip(*[libKeldysh.FieldToIntensity(Efield1), Ne_final]))
    lists = sorted(zip(*[libKeldysh.FieldToIntensity(Efield1), Ne_final, pulses1]))
    Intensity_sorted, Ne_final_sorted, pulses1_sorted = list(zip(*lists))
    plt.plot(Intensity_sorted, Ne_final_sorted/(pulses1), '-') #, label=r'$N_{exc}/\tau$') #, I='+str(libKeldysh.FieldToIntensity(Efield1_value))+r' W/m$^2$')
  
  ax12.set_yscale('log')
  ax12.set_ylabel(r'$w_{PI}$ (m$^{-3}$.s$^{-1}$)')

  #ax2 = plt.subplot(plotB)
  #plt.title('Virtual excitation factor')
  
  #plt.plot(libKeldysh.FieldToIntensity(Efield1), irreversibility, 'ow', label='Virtual exc.')
  
  #plt.xlabel(r'Peak intensity (W/m$^2$)')
  #plt.ylabel('Virtual exc. factor')
  #plt.ylim((0,1))
  #ax2.set_xscale('log')
  ##adding the Keldysh adiabatic factor
  #ax22 = ax2.twinx()
  
  #plt.plot(libKeldysh.FieldToIntensity(Efield1), libKeldysh.gammaKeldysh(2.56*e, 0.228, Efield1, wavelengths1), 'or', label='Adiab. fact.') #BUG: 2.56 must not be specified by hand here! 

  #plt.ylabel('Adiabadicity factor')
  ##plt.ylim((0,1))
  #ax22.set_yscale('log')

  #ax2.set_yscale("log")
  #ax2.legend(loc=2, borderaxespad=0.)
  #ax1.legend(bbox_to_anchor=(1.15, 1), loc=2, borderaxespad=0.)
  #ax22.legend(loc=1, borderaxespad=0.)

  # Saving the plot
  ExportName = "EffectOfIntensityPerCEP"
  plt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.eps')
  plt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.png')
  if(ShowFinalPlot==1):
    plt.show()

  del ax1, PlottingSinglePulses, wavelength_uniq
#}}}

## Nexc ( CEP )  
def plotEffectOfCEP(PlottingTable, pulseWavelength=800e-9, pulseDuration=20e-15, FieldAmplitude=1E9, PulseShape=3, pulseRiseTime=0, ShowFinalPlot=0, Egap_SI=2.56*e): #{{{
  print("#====== Plot: SINGLE PULSE, EFFECT OF THE CEP, SAME WAVELENGTHS, SAME POLARIZATIONS: final density as function of field ======")
  #Filter: double pulses, same wavelengths, same polarization

  # Simple pulse, single color
  PlottingSinglePulses = PlottingTable #get the full database
  PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, libPlotting.IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
  PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, libPlotting.IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0., libPlotting.IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, pulseWavelength, libPlotting.IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts. #TODO: this field should be making a loop in available wavelengths given in wavelengths1_uniq.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0, libPlotting.IndexTable['FreezeHXC']) #keeps only non-frozen XC.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, PulseShape, libPlotting.IndexTable['PulseType']) #keeps only non-frozen XC.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, pulseDuration, libPlotting.IndexTable['pulses1']) #removes fields where tau != pulseDuration
  PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, FieldAmplitude*1.01, libPlotting.IndexTable['Efield1'])
  PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, FieldAmplitude*0.99, libPlotting.IndexTable['Efield1'])
  if(pulseRiseTime != 0): 
    PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, pulseRiseTime, libPlotting.IndexTable['RiseTime']) #removes fields where tau != pulseDuration

  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = libPlotting.CountTheNumberOfLaserParameters(PlottingSinglePulses)
  
  plt.figure(figsize=(SizeX,SizeY))
  ax1 = plt.subplot(111) #111 or plotA if using plotB
  plt.title(r'$n_{\mathrm{exc}}(I)$, $\tau=$'+str(int(pulseDuration*1E15))+r' fs, $E_{peak}=$'+str(round(FieldAmplitude*1E-9))+' V/nm.')
  # Number of sweeped parameters: 3 (Efield, Delay, CrossPolarization)
  #plt.scatter(libKeldysh.FieldToIntensity(Efield1), Ne_final, c=Delay, s=(30.+CrossPolarization*50.), label=r'$N_{exc}$[Field]')
  plot1 = [] #np.empty((0,0))
  #for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
  #PlottingSinglePulses_CEP = libDatabase.FilterDatabase(PlottingSinglePulses, CEP1_value, libPlotting.IndexTable['CEP1']) #select the local CEP
  
  ## prints the current state of selection and stores the values of sweeped parameters
  runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingSinglePulses)

  
  #plot_field, = ax1.plot(libKeldysh.FieldToIntensity(Efield1), Ne_final, '-', label=r'$N_{exc}^{final}$, $\phi_{CE}$='+str(round(CEP1_value*180./pi)))
  #lists = sorted(itertools.izip(*[libKeldysh.FieldToIntensity(Efield1), Ne_final]))
  lists = sorted(zip(*[CEP1, Ne_final, pulses1]))
  CEP1_sorted, Ne_final_sorted, pulses1_sorted = list(zip(*lists))
  #print Header+"Ready to plot: "
  #print "** Info: Sorting data..."
  #print Efield1_sorted, Ne_final_sorted
  plot_field, = ax1.plot(np.multiply(CEP1_sorted,180./pi), Ne_final_sorted, 'o', label=r'$N_{exc}^{final}$')
  #print Header+"plotEffectOfIntensity_CEP(): adding CEP "+str(round(CEP1_value*180./pi))+" to the the plot."
  logger.info(plot_field)
  logger.debug(Header+"{Fore.RED}** Debug info: "+str(np.shape(plot_field)))
  plot1 += [plot_field] # does cumulate each CEP now. Type: Line2D, not np.array. 
    
    #ax1.legend(loc=4, borderaxespad=0.) #bbox_to_anchor=(1.05, 1), 
  labels1 = [l.get_label() for l in plot1]
  ax1.legend(plot1, labels1, loc='best', fancybox=True, framealpha=0)
  #plt.xlabel(r'Peak intensity (W/m$^2$)')
  plt.xlabel(r'Carrier envelope phase $\phi$ (deg)')
  plt.ylabel(r'Final density (m$^{-3}$)')
  #  plt.ylim((1E16,4.*5e28))
  #ax1.set_xscale('log')
  #ax1.set_yscale('log')
  #plot1_list = plot1.tolist()
    #print Efield1, CEP1_value
    
  #Adding axis of gammaKeldysh, but for all CEPs !!
  #runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, ElectronEnergy_max, Energy_max, W_pi, PulseType, RiseTime = libPlotting.ExtractData(PlottingSinglePulses)
  
  ax12 = ax1.twinx()
  #for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
  #PlottingSinglePulses_CEP = libDatabase.FilterDatabase(PlottingSinglePulses, CEP1_value, IndexTable['CEP1']) #select the local CEP
  
  ## prints the current state of selection and stores the values of sweeped parameters
  runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingTable)
  
  #plot_field, = ax1.plot(libKeldysh.FieldToIntensity(Efield1), Ne_final, '-', label=r'$N_{exc}^{final}$, $\phi_{CE}$='+str(round(CEP1_value*180./pi)))
  #lists = sorted(itertools.izip(*[libKeldysh.FieldToIntensity(Efield1), Ne_final]))
  lists = sorted(zip(*[CEP1, Ne_final, pulses1]))
  CEP1_sorted, Ne_final_sorted, pulses1_sorted = list(zip(*lists))
  plt.plot(np.multiply(CEP1_sorted,180./pi), np.divide(Ne_final_sorted, pulses1_sorted), '-') #, label=r'$N_{exc}/\tau$') #, I='+str(libKeldysh.FieldToIntensity(Efield1_value))+r' W/m$^2$')
  
  #ax12.set_yscale('log')
  ax12.set_ylabel(r'$w_{PI}$ (m$^{-3}$.s$^{-1}$)')

  #ax2 = plt.subplot(plotB)
  #plt.title('Virtual excitation factor')
  
  #plt.plot(libKeldysh.FieldToIntensity(Efield1), irreversibility, 'ow', label='Virtual exc.')
  
  #plt.xlabel(r'Peak intensity (W/m$^2$)')
  #plt.ylabel('Virtual exc. factor')
  #plt.ylim((0,1))
  #ax2.set_xscale('log')
  ##adding the Keldysh adiabatic factor
  #ax22 = ax2.twinx()
  
  #plt.plot(libKeldysh.FieldToIntensity(Efield1), libKeldysh.gammaKeldysh(2.56*e, 0.228, Efield1, wavelengths1), 'or', label='Adiab. fact.') #BUG: 2.56 must not be specified by hand here! 

  #plt.ylabel('Adiabadicity factor')
  ##plt.ylim((0,1))
  #ax22.set_yscale('log')

  #ax2.set_yscale("log")
  #ax2.legend(loc=2, borderaxespad=0.)
  #ax1.legend(bbox_to_anchor=(1.15, 1), loc=2, borderaxespad=0.)
  #ax22.legend(loc=1, borderaxespad=0.)

  # Saving the plot
  plt.tight_layout()
  ExportName = "EffectOfCEP"
  plt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.eps')
  plt.savefig('Metaplot-SinglePulse-SingleColor-'+ExportName+'.png')
  if(ShowFinalPlot==1):
    plt.show()

  del ax1, PlottingSinglePulses, wavelength_uniq
#}}}

## (Final density | Virtual excitation rate) (Intensity) with comparing FreezeHXC if available.
def plotEffectOfFrozenXC(PlottingTable): #{{{
  global meff
  print("#====== Plot: SINGLE PULSE, EFFECT OF THE FROZEN HXC, SAME WAVELENGTHS, SAME POLARIZATIONS: final density as function of field ======")
  #Filter: double pulses, same wavelengths, same polarization

  # Simple pulse, single color
  PlottingSinglePulses = PlottingTable #get the full database
  PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, libPlotting.IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
  PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, libPlotting.IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0., libPlotting.IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 800e-9, libPlotting.IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts. #TODO: this field should be making a loop in available wavelengths given in wavelengths1_uniq.
  #PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0, libPlotting.IndexTable['CEP1']) #keeps only CEP=0.


  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = libPlotting.CountTheNumberOfLaserParameters(PlottingSinglePulses)

  plt.figure()
  ax1 = plt.subplot(plotA)
  plt.title('Final density as function of peak intensity')
  # Number of sweeped parameters: 3 (Efield, Delay, CrossPolarization)
  #plt.scatter(libKeldysh.FieldToIntensity(Efield1), Ne_final, c=Delay, s=(30.+CrossPolarization*50.), label=r'$N_{exc}$[Field]')
  #for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
  #PlottingSinglePulses_CEP = libDatabase.FilterDatabase(PlottingSinglePulses, CEP1_value, libPlotting.IndexTable['CEP1']) #select the local CEP
  
  # Filtering the FreezeHXC states
  PlottingSinglePulses_CEP_NoFreezeHXC = libDatabase.FilterDatabase(PlottingSinglePulses, 0, libPlotting.IndexTable['FreezeHXC']) #select the non-frozen HXC
  PlottingSinglePulses_CEP_FreezeHXC = libDatabase.FilterDatabase(PlottingSinglePulses, 1, libPlotting.IndexTable['FreezeHXC']) #select the frozen HXC
  ## prints the current state of selection and stores the values of sweeped parameters



  runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingTable)


  plt.plot(libKeldysh.FieldToIntensity(Efield1), Ne_final, 'ow', label=r'$N_{exc}^{final}$') #, CEP='+str(round(CEP1_value*180./pi)))
  
  runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingTable)
  
  plt.plot(libKeldysh.FieldToIntensity(Efield1), Ne_final, 'ok', label=r'$N_{exc}^{final}$, FrozenHXC')#, CEP='+str(round(CEP1_value*180./pi)))
  
  ax1.legend(loc='best', fancybox=True, framealpha=0) #bbox_to_anchor=(1.05, 1), 
  #plt.xlabel(r'Peak intensity (W/m$^2$)')
  plt.ylabel(r'Final density (m$^{-3}$)')
  #  plt.ylim((1E16,4.*5e28))
  ax1.set_xscale('log')
  ax1.set_yscale('log')

  #Adding axis of gammaKeldysh
  ax12 = ax1.twinx()
  
  runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max,Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingSinglePulses_CEP_NoFreezeHXC)

  
  plt.plot(libKeldysh.FieldToIntensity(Efield1), Ne_final/pulses1, 'ow') #, label=r'$N_{exc}/\tau$') #, I='+str(libKeldysh.FieldToIntensity(Efield1_value))+r' W/m$^2$')
  



  runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingTable)

  
  plt.plot(libKeldysh.FieldToIntensity(Efield1), Ne_final/pulses1, 'ok') #, label=r'$N_{exc}/\tau$, FrozenHXC') #, I='+str(libKeldysh.FieldToIntensity(Efield1_value))+r'
  
  ax12.set_yscale('log')
  ax12.set_ylabel(r'$w_{PI}$ (m$^{-3}$.s$^{-1}$)')

  ax2 = plt.subplot(plotB)
  plt.title(r'Virtual excitation $\zeta$')
  



  runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingTable)

  
  plt.plot(libKeldysh.FieldToIntensity(Efield1), irreversibility, 'ow', label='Virtual exc.')
  



  runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingTable)


  plt.plot(libKeldysh.FieldToIntensity(Efield1), irreversibility, 'ok', label='Virtual exc., FreezeHXC')
  
  plt.xlabel(r'Peak intensity (W/cm$^2$)')
  plt.ylabel(r'Virtual exc. $\zeta$')
  plt.ylim((0,1))
  ax2.set_xscale('log')
  #adding the Keldysh adiabatic factor
  ax22 = ax2.twinx()
  



  runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingSinglePulses_CEP_FreezeHXC)

  
  plt.plot(libKeldysh.FieldToIntensity(Efield1), libKeldysh.gammaKeldysh(Egap_SI, meff, Efield1, wavelengths1), 'or', label=r'Adiabaticity $\gamma$') #BUG: 2.56 must not be specified by hand here! 
  



  runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingSinglePulses_CEP_FreezeHXC)

  
  plt.plot(libKeldysh.FieldToIntensity(Efield1), libKeldysh.gammaKeldysh(Egap_SI, meff, Efield1, wavelengths1), 'or') #, label='Adiabaticity') #BUG: 2.56 must not be specified by hand here! 
  
  plt.ylabel(r'Adiabadicity $\gamma$')
  #plt.ylim((0,1))
  ax22.set_yscale('log')

  #print Efield1, CEP1_value
  # ax2.set_yscale("log")
  #ax2.legend(loc=2, borderaxespad=0.)
  # ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
  #ax22.legend(loc=1, borderaxespad=0.)
  plt.tight_layout()
  # Saving the plot
  plt.savefig('Metaplot-SinglePulse-SingleColor-EffectOfFrozenHXC.eps')
  plt.savefig('Metaplot-SinglePulse-SingleColor-EffectOfFrozenHXC.png')
  #if(ShowFinalPlot==1):
    #plt.show()

  del ax1, PlottingSinglePulses, wavelength_uniq, wavelengths1
#}}}

