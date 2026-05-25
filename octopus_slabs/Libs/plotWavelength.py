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

## @module plotWavelength
# Plots showing wavelength-dependent observables (Nexc, Wpi, Energy)
# for single- and two-colour irradiations.
# Functions: plotEffectOfWavelength_Tolerance,
#            Bicolor_plotEffectOfWavelength_Tolerance,
#            Bicolor_plotEffectOfCEP,
#            plotEffectOfWavelength_2dMap

from octopus_slabs.Libs._plotCommon import *  # noqa: F401,F403

## Plots (Nexc | Wpi | Energy) (wavelength) (1D plot) by introducing <FieldToTolerance> percents tolerance on the electric field. 
# This enables to analyze the role of laser wavelength on the final excited electron density, excitation rate and laser-induced energy increase. Provides predictions on the Stark effect. 
# Wpi is compared to analytical Keldysh theory. 
# Energy is compared with analytical Drude model. 
def plotEffectOfWavelength_Tolerance(PlottingTable, FieldTolerance=0.05, pulseDuration=20e-15, CEP=0., PulseType=3,
                                     RiseTimeRequested=2, ShowFinalPlot=0, NumberOfAtoms=2., Egap_SI=2.56*e,
                                     Egap_SI_indirect=0.51*e, wavelengths1_min=200e-9, wavelengths1_max=3300e-9,
                                     ShowDensities=False, PlotIrreversibility=False, PlotFreeCarrierAbsorption=True,
                                     EnergyAxis=True, Functional="LDA"): #{{{
  PlotStarkShift  = False
  ShowInset       = False
  PlotValleysInset= True #True: plot valleys, false: plot hills (max of Wpi)
  PlotSaturation  = False
  AccountForOpticalIndex=False #True: account for linear optical spectrum. False: does not account. 
  #WARNING: TDLDA optical spectrum from Reining and Rubio is not reaching our frequencies of interest!

  global n0
  ContinuousKeldysh_Nb_Fields     =100  #100: not really sufficient, 500: Good for publication.
  ContinuousKeldysh_Nb_Wavelengths=500  #500: not sufficient.        3000: Good for publication.
  Use128bitsFloats=False; 
  OpticalDatabase="Palik" #"TDLDA" (down to 2 eV, then extrapolation, very imprecise!! | "Palik" (super wide)
  InterpolationMethod="nearest" #"UnivariateSpline", "cubic", "linear", 'nearest'
  DrawOrders=6 #only visual, nothing about precision
  
  meff_Waroquiers_LDA = 0.2226 #effective mass of electron
  #[Laflamme Janssen, J.; Gillet, Y.; Poncé, S.; Martin, A.; Torrent, M. & Gonze, X. Precise effective masses from density functional perturbation theory Physical Review B, 2016, 93, 205147]
  #meff = 0.14 #0.14: fits for Si in MPI regime at 800 nm.
  #Field dependent fitting of effective mass
#  meff_LowField = meff_Waroquiers_LDA #0.14
  meff_HighField = meff_Waroquiers_LDA #0.14: fits for Si in MPI regime at 800 nm.
#  meff_VeryHighField = meff_Waroquiers_LDA #0.5
  meff = meff_Waroquiers_LDA
  meff_Tunnel = meff_Waroquiers_LDA #fits well if using INDIRECT band gap for 800 nm
  print("#====== Plot: SINGLE PULSE, EFFECT OF THE WAVELENGTH, SAME CEP, SAME POLARIZATIONS: final density as function of field ======")
  #Filter: double pulses, same wavelengths, same polarization
  
  LineWidthCoeff_abs = 0.5 #controls absolute thickness of all curves
  LineWidthCoeff_rel = 1.0 #controls thickness of Keldysh-Gruzdev curves relatively to TDDFT
  antialias=True
  
  #FieldOfInterest = [5e8, 1e9, 2e9, 3E9, 4E9, 5e9, 1e10]
  FieldOfInterest = [5e8, 7.5E8, 1e9, 1.5e9, 2e9, 3e9, 4E9, 5e9] #, 1e10]
  NumberOfColors=len(FieldOfInterest)
  ColorList=[plt.cm.seismic(i) for i in np.linspace(0, 1, NumberOfColors)] 
  # Not bad! RdYlBu, coolwarm, seismic
  # Too soft: Copper, Spectral
  # Not perfect. RdBu. RdYlGn.
  # Too dynamic: bwr
  # To test: , brg, rainbow
  #Colors={5e8: 'c', 7.5e8: 'r', 1e9: 'b', 1.5e9: 'g', 2e9: 'y', 3e9: 'm'} #, 1e10: 'k'}
  Colors={FieldOfInterest[i]:ColorList[i] for i in range(0,NumberOfColors)}
  LineWidths={FieldOfInterest[i]:np.arange(1.5,4,NumberOfColors) for i in range(0,NumberOfColors)}
  #EffectiveMass={1e9: meff_LowField, 2e9: meff_LowField, 5e9: meff_HighField, 1e10: meff_VeryHighField}
  EffectiveMass={FieldOfInterest[i]:meff_HighField for i in range(0,NumberOfColors)}
  # Simple pulse, single color
  


  ### TODO: MOVE TO libPlotting.py to access it for plotTimeKeldysh and plotTimeHHG.py. 
  def ImportPalikDatabase_epsilon_twofiles_eV(wavelengths, DataFile='Database/Si-TDLDA'): #{{{
                                                                                                  
    Header="** ImportPalikDatabase_epsilon_twofiles: "
    import os
    DataFile=os.environ["spp_extended_theory"]+"/spp_extended_theory/SimpleSPP/Database/Si-TDLDA"
    nfile = DataFile+"-epsR.csv" #TODO: rename nfile to ReEpsFile
    kfile = DataFile+"-epsC.csv" #TODO: rename kfile to ImEpsFile
    
    logger.debug(str(nfile)+", "+str(kfile))
    
    #print("** Info: opening "+nfile+" and "+kfile+".")
    narray = np.loadtxt(nfile, delimiter="\t", skiprows=2)
    karray = np.loadtxt(kfile, delimiter="\t", skiprows=2)
    
    unit = 1E0 #Unit of the wavelength found in databases <nfile> and <kfile>. 
    
    # import wavelength, epsilonR, epsilonC
    energies1 = narray[:,0]*unit #/!\ its an energy!
    energies2 = karray[:,0]*unit #/!\ its an energy!
    n = narray[:,1]
    kk = karray[:,1]
    #numrows = 10000 #Number of rows to interpolate the data on. 
    base = 10
    # interpolate on new wavelength mesh using 1st order
    order=1 #linear interpolation
    #print Header+"** Info: Generating the new wavelength mesh: ("+str(np.amin(energies2))+", "+str(np.amax(energies2))+")"
    #wavelengths = np.logspace(np.amin(np.log10(energies2)), np.amax(np.log10(energies2)), num=numrows, base=base, endpoint = True)
    
    #print Header+"New wavelength mesh has "+str(numrows)+" rows."
    #print wavelengths
    
    #print(Header+"** Info: definition of interpolation functions...")
    energy1=np.sort(energies1)
    energy2=np.sort(energies2)
    del energies1, energies2 
    if(InterpolationMethod=='UnivariateSpline'):
      fni = InterpolatedUnivariateSpline(energy1, n, k=order)
      fki = InterpolatedUnivariateSpline(energy2, kk, k=order)
    elif(InterpolationMethod=='cubic'):
      fni = interp1d(energy1, n, kind='cubic', fill_value='extrapolate')
      fki = interp1d(energy2, kk, kind='cubic', fill_value='extrapolate')
    elif(InterpolationMethod=='linear'):
      fni = interp1d(energy1, n,  kind='linear', fill_value='extrapolate')
      fki = interp1d(energy2, kk, kind='linear', fill_value='extrapolate')
    else: 
      fni = interp1d(energy1, n, kind='nearest', fill_value='extrapolate')
      fki = interp1d(energy2, kk, kind='nearest', fill_value='extrapolate')
    ##Interpolation for one optical constant
    #try:
      #ni = fni(wavelength); ki = fki(wavelength)
      ##print wavelength, ni, ki
      #epsilon = (ni+1j*ki) #NOTE: we are picking up the epsRe, and epsIm directly here
      #print "" 
      #print "Interpolated permittivity at "+str(wavelength*1E9)+" nm = "+str(epsilon)
    #except: 
      #print "Interpolation for "+str(wavelength*1E9)+" nm failed."
	    
	    
    # We need to express the required energy mesh into a wavelength
    energies = hbar*np.divide(2.*np.pi * c, wavelengths) / e #eV
    
    ## WARNING: Names are misleading here: we actually work with dielectric permittivities! 
    # defining the new epsR and epsC on a common mesh
    ni = fni(energies) #eV
    ki = fki(energies)
    epsilons = np.add(ni, np.multiply(1.j, ki)) 
    return epsilons
  #}}}
  
  
  
  # TODO: Filtering function
  PlottingSinglePulses = PlottingTable #get the full database
  PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, libPlotting.IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
  PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, libPlotting.IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0., libPlotting.IndexTable['CrossPolarization']) #removes fields where CrossPolarization is different from 0.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, PulseType, libPlotting.IndexTable['PulseType']) #keeps only pulses of type StH (3)
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, libPlotting.IndexFunctional[Functional], libPlotting.IndexTable['Functional']) #removes runs where functional is not LDA
  #PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, RiseTime, libPlotting.IndexTable['RiseTime']) #keeps only pulses of type STH RiseTime
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, pulseDuration, libPlotting.IndexTable['pulses1']) #keeps only pulses of selected duration in argument

  if(RiseTimeRequested == 2): #select by RiseTime only if it is defined in parameters
      PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, RiseTimeRequested, libPlotting.IndexTable['RiseTime']) #but this removes some pulses that give 0 althought they do have 2 cycle risetime.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, CEP, libPlotting.IndexTable['CEP1']) #keeps only CEP = 0.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0., libPlotting.IndexTable['FieldInduced']) #remove fields computed with InducedFields enabled.

  PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, wavelengths1_max, libPlotting.IndexTable['wavelengths1']) #we will try to add data from 3200 nm.

  ## TODO: check if CEP could count at smaller wavelength. 

  logger.info(Header+"Available data: ")
  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = libPlotting.CountTheNumberOfLaserParameters(PlottingSinglePulses)

  #if(len(CEP1_uniq) == 1): #{{{ #single color simulations
  # Complexity level: 2
  
  if(output_format=='publication'):
    plt.figure(figsize=(6,2*5)) #new figure
    plotA=211; plotB=212 #vertical figures (adequate for PRL)
  else: 
    plt.figure(figsize=(2*6,6)) #new figure
    plotA=121; plotB=122 #horizontal figures (adequate for presentation)
  
  ax1 = plt.subplot(plotB)
  #ax1.set_prop_cycle('color',Colors)
  
  # Create inset in ax1 region. 
  if(ShowInset):
    ax11 = inset_axes(ax1, width="37%", height="30%", loc="lower right", borderpad=1.3)
  #ax11 = inset_axes([0.5, 0.05, 0.95, 0.5])
  #ax11.set_xmargin(0.5)
  #ax11.set_ymargin(2)
  # Turning on labels for the inset
    ax11.tick_params(top=True, left=True, bottom=False, labelleft=True, labeltop=True, labelbottom=False)
  #ax11.tick_params(axis='both', top=True, left=True)
  #wpi_peaks_wavelengths  =np.zeros((len(FieldOfInterest),20)); 
  #wpi_valleys_wavelengths=np.zeros((len(FieldOfInterest),20)); 
  #wpi_peaks_Efields      =np.zeros((len(FieldOfInterest),20)); 
  #wpi_valleys_Efields    =np.zeros((len(FieldOfInterest),20)); 
  #wpi_peaks_values       =np.zeros((len(FieldOfInterest),20)); 
  #wpi_valleys_values     =np.zeros((len(FieldOfInterest),20));
  
  Peaks=[]; Valleys=[]
  Keldysh_Peaks=[]; Keldysh_Valleys=[]
  
  plot1 = []; plot11=[] #np.empty((0,0))

  ## Returns position and value at the peaks and valley.
  def FindMinAndMax(signal, width):  # {{{
      Header = "plotFinalQuantities: FindMinAndMax():"
      logger.info(Header)
      # print np.shape(signal), np.shape(width)
      width = np.asfarray(width)
      # print width
      if (EnablePeakFinder):
          peaks, prop_peaks = find_peaks(signal)  # , width)
          valleys, prop_valleys = find_peaks(-signal)
      else:
          peaks = [];
          prop_peaks = []
          valleys = [];
          prop_valleys = []
      # print peaks, signal[peaks], width[peaks]

      # print "Min: "+str(valleys)
      # print "Max: "+str(peaks)
      return width[peaks], signal[peaks], width[valleys], signal[valleys]

  #}}}

  for Efield1_value in Efield1_uniq: #for each value of Efield1, we make a new search, and plot one new line in the chart #{{{
    #PlottingSinglePulses_Efield = libDatabase.FilterDatabase(PlottingSinglePulses, Efield1_value, IndexTable['Efield1']) #select the local Efield




    PlottingSinglePulses_Efield = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, Efield1_value*(1.+FieldTolerance), libPlotting.IndexTable['Efield1'])

    PlottingSinglePulses_Efield = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses_Efield, Efield1_value*(1.-FieldTolerance), libPlotting.IndexTable['Efield1'])
    #print Header+"Dataset for the field: "+str(Efield1_value*1e-9)+' V/nm'
    #print PlottingSinglePulses_Efield
    # TODO: as fields are all slighly different, we should project on the neareast neighboring value, up to 1E8 V/m. 
    
    ## prints the current state of selection and stores the values of sweeped parameters



    runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingSinglePulses_Efield)

    
    # Replaces wavelengths with photon energy
    if(EnergyAxis):  #very dirty
       wavelengths1_Show = np.divide(2*np.pi*c*hbar/e, wavelengths1) #eV #E=hbar omega = 2 pi c hbar / wavelength
    else:
       wavelengths1_Show = wavelengths1
    
    ## Sorting data before plotting. 
    # If using Energies, then we sort by energy. 
    # If using wavelengths, sorting is made by wavelengths. 
    lists = sorted(zip(*[wavelengths1_Show, Ne_final, pulses1, ElectronEnergy_max]))
    wavelengths1_sorted, Ne_final_sorted, pulses1_sorted, Energy_max_sorted = list(zip(*lists))
    
    if(not EnergyAxis):
        wavelengths1_sorted_show = np.multiply(1e9, wavelengths1_sorted)
    else: #if our wavelength is actually an energy
        wavelengths1_sorted_show = wavelengths1_sorted
    wpi_sorted = np.divide(Ne_final_sorted, np.multiply(1e0, pulses1_sorted))
    
    ## Returns position and value at the peaks and valley. 
    def FindMinAndMax(signal, width): #{{{
        Header="plotFinalQuantities: FindMinAndMax():"
        logger.info(Header)
        #print np.shape(signal), np.shape(width)
        width=np.asfarray(width)
        #print width
        if(EnablePeakFinder): 
          peaks, prop_peaks = find_peaks(signal) #, width)
          valleys, prop_valleys= find_peaks(-signal)
        else: 
          peaks=[]; prop_peaks=[]
          valleys=[]; prop_valleys=[]
        #print peaks, signal[peaks], width[peaks]
 
        #print "Min: "+str(valleys)
        #print "Max: "+str(peaks)
        return width[peaks], signal[peaks], width[valleys], signal[valleys]
    #}}}   
    #array_id=0
    for Efield1_compare in FieldOfInterest: #{{{ #WARNING: Efield1_value should be considered randomized
    #TODO: PARALLELIZE THIS LOOP. It can be VERY long. 
        if(abs(Efield1_value - Efield1_compare)/Efield1_value < FieldTolerance): #{{{
            wavelengths_peaks, wpi_peaks, wavelengths_valleys, wpi_valleys=FindMinAndMax(wpi_sorted, wavelengths1_sorted)
            logger.info("W_PI peaks are at: "+str(wavelengths_peaks))
            logger.info("W_PI hills are at: "+str(wavelengths_valleys))
            NumberOfPeaks   = np.shape(wavelengths_peaks)
            NumberOfValleys = np.shape(wavelengths_valleys)
            logger.info("Adding extrema to a list for plotting...")
            for element in wavelengths_peaks:
                Peaks.append(  [Efield1_compare, element]) #one line per datapoint
            for element in wavelengths_valleys:
                Valleys.append([Efield1_compare, element]) #one line per datapoint
            #wpi_peaks_wavelengths  [array_id]=wavelengths_peaks  
            #wpi_peaks_values       [array_id]=wpi_peaks          
            #wpi_peaks_Efields      [array_id]=Efield1_compare    
            #wpi_valleys_wavelengths[array_id]=wavelengths_valleys
            #wpi_valleys_values     [array_id]=wpi_valleys        
            #wpi_valleys_Efields    [array_id]=Efield1_compare    
            
            #array_id+=1 #we need something to localize the values
        #### PLOTTING
            #n0 = 4.*Density
            #if(not EnergyAxis):
            plot_field, = ax1.plot(wavelengths1_sorted_show, wpi_sorted,'-^', antialiased=antialias, color=Colors[Efield1_compare], linewidth=LineWidthCoeff_abs*LineWidths[Efield1_compare], label=r'$\widetilde{w}_{\mathrm{PI}}$, I='+str(libKeldysh.FieldToIntensity(Efield1_value))+r' W/m$^2$')
            
            #plot_field, = ax1.plot(wavelengths1_sorted_show, wpi_sorted,'-^', linewidth=LineWidthCoeff_abs*LineWidths[Efield1_compare], label=r'$\widetilde{w}_{\mathrm{PI}}$, I='+str(libKeldysh.FieldToIntensity(Efield1_compare))+r' W/m$^2$')
            plot1 += [plot_field]
            
            #}}} #end if
            
    #}}} #end for
    del PlottingSinglePulses_Efield
  #}}} #end for
        
  Peaks=np.array(Peaks).T; Valleys=np.array(Valleys).T
  logger.info(str(np.shape(Peaks))+", "+str(np.shape(Valleys)))
  ##Peaks_wavelengths=Peaks[0]
  if(ShowInset):
    if(PlotValleysInset):
        ax11.scatter(np.multiply(1E-9,Valleys[0,:]), Valleys[1,:], s=10, c="red", marker='^', label="Valleys TDDFT")
    else: 
        ax11.scatter(np.multiply(1E-9,Peaks[0,:]), Peaks[1,:], s=10, c="red", marker='^', label="Peaks TDDFT")
  #plot11+=[plot_field]
  #ax11.plot(wpi_peaks_Efields,   wpi_peaks_wavelengths  )
  
  #exit()
  # We separate TDDFT and Keldysh predictions now. 
  ax2 = plt.subplot(plotA)
  # Create inset in ax2 region. 
  ax22 = inset_axes(ax2, width="35%", height="35%", loc="lower right", borderpad=1.5)
  
  # Turning on labels for the inset
  ax22.tick_params(labelleft=True, labelbottom=True)
  
  ## Preparing the data for Keldysh computations
  order = 50  #Keldysh series order
  
  if(PlotStarkShift): #{{{
    # Can we plot Wpi_ResonantGap(EgapEff) ? 
    IntegerLevels=np.arange(1,6,1) #generates 1,2,3,...,10
    for IntegerLevel in IntegerLevels:
        EgapResonant=Egap_SI/IntegerLevel #these are my levels resonant with initial gap
        EgapResonant_Wavelength=h*c/EgapResonant #gives resonant wavelength with field
        
        epsilon_Si0_Stark = ImportPalikDatabase_epsilon_twofiles_eV(EgapResonant_Wavelength, "Database/Si-TDLDA")
        OpticalIndex_Stark = np.real(np.sqrt(epsilon_Si0_Stark)) #maybe problem on the sqrt(real(sqrt(eps)))?
        
        # We compute Wpi_Gruzdev at this multiply resonant wavelength
        Wpi_ResonantGap = libKeldysh.GenerateKeldyshDatabase(Egap_SI, meff_Waroquiers_LDA, EgapResonant_Wavelength, np.multiply(FieldOfInterest, np.sqrt(OpticalIndex_Stark)), order, 1.) #NOTE that OpticalIndex*FieldOfInterest was ignored for now. 
        gammaResonant = libKeldysh.gammaKeldysh(Egap_SI, meff_Waroquiers_LDA, np.multiply(FieldOfInterest, np.sqrt(OpticalIndex_Stark)), EgapResonant_Wavelength, 1.) #valid for scalar data
        k1 = libKeldysh.Keldysh1phi(gammaResonant); k2 = libKeldysh.Keldysh2theta(gammaResonant) #valid
        EgapEff = libKeldysh.EffectiveGap(Egap_SI, k1, k2)
        del k1, k2

        ## We re-inject the new resonant wavelengths into Keldysh, to compute the "real Wpi"
        #Wpi_ResonantGap = libKeldysh.GenerateKeldyshDatabase(Egap_SI, meff_Waroquiers_LDA, h*c/EgapEff, np.multiply(FieldOfInterest, np.sqrt(OpticalIndex_Stark)), order, 1.) #NOTE that
        
        plot_field, = ax1.plot((EgapEff/e/IntegerLevel), Wpi_ResonantGap, '--k', antialiased=antialias, linewidth=1)
        plot1 += [plot_field]
        plot_field, = ax2.plot((EgapEff/e/IntegerLevel), Wpi_ResonantGap, '--k', antialiased=antialias, linewidth=1)
        plot1 += [plot_field]
        del EgapEff, Wpi_ResonantGap
  #}}}
  

  #wavelengths1_min = 200e-9; wavelengths1_max = 3300e-9 #2100e-9
  #wavelengths1_Keldysh_log = np.linspace(np.log10(np.float128(wavelengths1_min)), np.log10(np.float128(wavelengths1_max)), ContinuousKeldysh_Nb_Wavelengths)
  #wavelengths1_Keldysh     = np.power(10., wavelengths1_Keldysh_log)
  wavelengths1_Keldysh = np.linspace(wavelengths1_min, wavelengths1_max, ContinuousKeldysh_Nb_Wavelengths)
  #wavelengths1_min = 1e-6; wavelengths1_max = 20e-6 #3300e-9 #2100e-9

  ### Capturing the dielectric permittivity in experimental or TDLDA database
  if(OpticalDatabase == "Palik"):
    epsilon_Si0 = libOpticalData.ImportPalikDatabase_epsilon_fromNK(wavelengths1_Keldysh) ## WARNING: replacing by LDA is DIFFICULT, since PRL Rubio is defined down to 2 eV. Nothing below!!
  elif(OpticalDatabase == "TDLDA"):
    logger.warning(Header+"** WARNING! TDLDA database is poorly defined in mid-IR !!")
    epsilon_Si0 = ImportPalikDatabase_epsilon_twofiles_eV(wavelengths1_Keldysh, "Database/Si-TDLDA")
  OpticalIndex = np.real(np.sqrt(epsilon_Si0))
  OpticalIndex_disabled = np.divide(OpticalIndex, OpticalIndex)
  
  if(not AccountForOpticalIndex):
      OpticalIndex = OpticalIndex_disabled 
    
  logger.info(Header+"Building comparison with Keldysh-Gruzdev...")
  #wavelengths1_min = np.min(wavelengths1); wavelengths1_max = np.max(wavelengths1)
  logger.info(Header+"** INFO TEST: Dielectric permittivity of Si at 800 nm")
  logger.info("[Palik] | [Reining et al PRL 2002]")
  logger.info("400 nm: "+str(libOpticalData.ImportPalikDatabase_epsilon_fromNK(400e-9))+"   | "+str(ImportPalikDatabase_epsilon_twofiles_eV(400e-9, "Database/Si-TDLDA")))
  logger.info("800 nm: "+str(libOpticalData.ImportPalikDatabase_epsilon_fromNK(800e-9))+"   | "+str(ImportPalikDatabase_epsilon_twofiles_eV(800e-9, "Database/Si-TDLDA")))
  logger.info("3200 nm: "+str(libOpticalData.ImportPalikDatabase_epsilon_fromNK(3200e-9))+" | "+str(ImportPalikDatabase_epsilon_twofiles_eV(3200e-9, "Database/Si-TDLDA")))
    
         
  for Efield1_compare in FieldOfInterest: #{{{ #INFO: good for Fig3b. 
    #NOTE: actually we can build comparison using Keldysh for much more data points. Only TDDFT limits us in number of points. And there is nothing about TDDFT in this loop! 
    
    ### I will need intermerdiate Keldysh coefficients
    #gamma = libKeldysh.gammaKeldysh(Egap, meff, PeakField, wavelength, RefractiveIndex.real) #valid for scalar data
    #k1 = Keldysh1phi(gamma); k2 = Keldysh2theta(gamma)
    #EgapEff=EffectiveGap(Egap, k1, k2)
    
    ## We construct Wpi(corresponding to the photon energy equal to the effective gap), and plot that! 
    
    #PeakWpi_Gruzdev_Stark = 
    
    PeakWpi_Gruzdev = libKeldysh.GenerateKeldyshDatabase(Egap_SI, EffectiveMass[Efield1_compare], wavelengths1_Keldysh, Efield1_compare*np.sqrt(OpticalIndex), order, 1., False) #NOTE: Refractive index is accounted in field, not in adiabadicity coefficient
    PeakWpi_Gruzdev_frozen = libKeldysh.GenerateKeldyshDatabase(Egap_SI, EffectiveMass[Efield1_compare], wavelengths1_Keldysh, Efield1_compare*np.sqrt(OpticalIndex), order, 1., True) #NOTE: Refractive index is accounted in field, not in adiabadicity coefficient
    Saturation = (n0 - (PeakWpi_Gruzdev * pulseDuration )) / n0 #saturation term
    PeakWpi_Gruzdev_sat = PeakWpi_Gruzdev * Saturation
    #PeakWpi_Gulley = GenerateKeldyshGulleyDatabase(Egap_SI, EffectiveMass[Efield1_compare], wavelengths1_Keldysh, Efield1_value, order, 1.)
    #PeakWpi_Gruzdev_Indirect = libKeldysh.GenerateKeldyshDatabase(Egap_SI_indirect, meff, wavelengths1_Keldysh, Efield1_value, order)
    
    PeakWpi_Tunneling_Direct = np.clip(libKeldysh.KeldyshTunnelingLimit(Egap_SI, meff_Tunnel, wavelengths1_Keldysh, Efield1_compare, False, 1.), 1E36, 1E46) #NOTE: physically, refractive index does not affect tunneling time
    
    #PeakWpi_Tunneling_Indirect = np.clip(libKeldysh.KeldyshTunnelingLimit(Egap_SI_indirect, meff_Tunnel, wavelengths1_Keldysh, Efield1_value), 1E36, 1E46)
    #PeakWpi_Gruzdev=
    #PeakWpi_Effective = np.power(np.power(PeakWpi_Gruzdev, -1) + np.power(PeakWpi_Tunneling_Direct,-1),-1)
    # Dirty: choosing energy or wavelengths. Sorting the results with this choice. 
    if(EnergyAxis):
        wavelengths1_Keldysh_Show = np.divide(2.*np.pi*c*hbar/e, wavelengths1_Keldysh) #eV
        wavelengths1_Keldysh_frozen_Show = np.divide(2.*np.pi*c*hbar/e, wavelengths1_Keldysh) #eV
    else:
        wavelengths1_Keldysh_Show = np.multiply(1E9, wavelengths1_Keldysh)
        wavelengths1_Keldysh_frozen_Show = np.multiply(1E9, wavelengths1_Keldysh)
    
    # resorting the results to plot as function of energy if necessary
    lists = sorted(zip(*[wavelengths1_Keldysh_Show, PeakWpi_Gruzdev, PeakWpi_Gruzdev_sat]))
    wavelengths1_Keldysh_Show, PeakWpi_Gruzdev_Show, PeakWpi_Gruzdev_sat_Show = list(zip(*lists)) 
    
    lists_frozen = sorted(zip(*[wavelengths1_Keldysh_frozen_Show, PeakWpi_Gruzdev_frozen]))
    wavelengths1_Keldysh_frozen_Show, PeakWpi_Gruzdev_frozen_Show = list(zip(*lists_frozen)) 
    
    plot_field, = ax2.plot(wavelengths1_Keldysh_Show, PeakWpi_Gruzdev_Show, '-', antialiased=antialias, color=Colors[Efield1_compare], linewidth=LineWidthCoeff_rel*LineWidthCoeff_abs*LineWidths[Efield1_compare], label=r'$\widetilde{w}_{\mathrm{PI}}^{\mathrm{KG}}$, I='+str(libKeldysh.FieldToIntensity(Efield1_compare))+r' W/cm$^2$')
    plot1 += [plot_field]
    
    ## Plotting Keldysh-Gruzdev theory in inset of Fig (bottom) (ax22). 
    plot_field, = ax22.plot(wavelengths1_Keldysh_frozen_Show, PeakWpi_Gruzdev_frozen_Show, '-', antialiased=antialias, color=Colors[Efield1_compare], linewidth=LineWidthCoeff_rel*LineWidthCoeff_abs*LineWidths[Efield1_compare], label=r'$\widetilde{w}_{\mathrm{PI}}^{\mathrm{KG}}$ frozen, I='+str(libKeldysh.FieldToIntensity(Efield1_compare))+r' W/cm$^2$')
    
    
    if(PlotSaturation): 
        ## Let's plot here the saturation effect
        plot_field, = ax2.plot(wavelengths1_Keldysh_Show, PeakWpi_Gruzdev_sat_Show, 
                    '-.', color=Colors[Efield1_compare], linewidth=LineWidthCoeff_rel*LineWidthCoeff_abs*LineWidths[Efield1_compare], label=r'$w_{\mathrm{PI}}\frac{n_0-n_{\mathrm{exc}}{n_0}$, I='+str(libKeldysh.FieldToIntensity(Efield1_compare))+r' W/m$^2$')
        plot1 += [plot_field]
                           #np.multiply(
            #PeakWpi_Gruzdev_Show, np.divide(np.add(n0, np.multiply(-1,PeakWpi_Gruzdev_Show*pulses1_sorted[0])), n0)),
    
    
    
    #plot_field, = ax1.plot(wavelengths1_Keldysh_Show, PeakWpi_Gulley, '-.'+Colors[Efield1_compare], label=r'$\widetilde{w}_{\mathrm{PI}}$ [Gulley], I='+str(libKeldysh.FieldToIntensity(Efield1_value))+r' W/m$^2$')
    #plot1 += [plot_field]
    
    #plot_field, = ax1.plot(wavelengths1_Keldysh_Show, PeakWpi_Tunneling_Direct, '.-'+Colors[Efield1_compare], label=r'$\widetilde{w}_{\mathrm{TI}}^{\mathrm{(d)}}$, I='+str(libKeldysh.FieldToIntensity(Efield1_value))+r' W/m$^2$')
    #plot1 += [plot_field]
    #plot_field, = ax1.plot(wavelengths1_Keldysh_Show, PeakWpi_Tunneling_Indirect, '-.'+Colors[Efield1_compare], label=r'$\widetilde{w}_{\mathrm{TI}}^{\mathrm{(ind)}}$, I='+str(libKeldysh.FieldToIntensity(Efield1_value))+r' W/m$^2$')
    #plot1 += [plot_field]
    #plot_field, = ax1.plot(wavelengths1_Keldysh_Show, PeakWpi_Effective, '-.'+Colors[Efield1_compare], label=r'$\widetilde{w}_{\mathrm{eff}}$, I='+str(libKeldysh.FieldToIntensity(Efield1_value))+r' W/m$^2$')
    #plot1 += [plot_field]
    
    ## SENSELESS. 
    #plot_field, = ax1.plot(wavelengths1_Keldysh_Show, PeakWpi_Gruzdev_Indirect, '--.'+Colors[Efield1_compare], label=r'$\widetilde{w}_{\mathrm{PI}}$ [Gruzdev], I='+str(libKeldysh.FieldToIntensity(Efield1_value))+r' W/m$^2$')
    #plot1 += [plot_field]

    #plot_field, = ax1.plot(wavelengths1_Keldysh_Show, PeakWpi_Gruzdev_Show, '-', linewidth=LineWidthCoeff_rel*LineWidthCoeff_abs*LineWidths[Efield1_compare], label=r'$\widetilde{w}_{\mathrm{PI}}$ [Gruzdev], I='+str(libKeldysh.FieldToIntensity(Efield1_value))+r' W/m$^2$')
    
    plot1 += [plot_field]
  #}}} #end for plotting Fig2b
  
  # Main figure will account for field, but the inset should not account for field [why?]
  #if(not AccountForOpticalIndex):
  #INFO: provides the Keldysh analysis part for Fig3a (inset).
  #WARNING: here we have disabled the contribution of optical index. Then one obtain the match for gap resonances with TDDFT. 
  # However, it seems W_PI amplitude is good only if we account for optical index. How the optical index should be treated in Keldysh theory or in TDDFT ?
  if(not AccountForOpticalIndex): 
      CoeffForDrawingValuesBeyond = 4 #we simply plot more values on grey lines, since TDDFT data have higher W_PI than Keldysh (when disregarding optical index correction)
  else: 
      CoeffForDrawingValuesBeyond=1 #optical index is sufficient to cover all the needed points in the plot
  ContinuousKeldysh_Fields_Log=np.linspace(8,np.log10(CoeffForDrawingValuesBeyond*np.max(FieldOfInterest)*np.sqrt(np.max(OpticalIndex))),ContinuousKeldysh_Nb_Fields)
  ContinuousKeldysh_Fields    =np.power(10., ContinuousKeldysh_Fields_Log)
  
  for Efield1_compare in ContinuousKeldysh_Fields:
    
    PeakWpi_Gruzdev_continuous = libKeldysh.GenerateKeldyshDatabase(Egap_SI, meff_Waroquiers_LDA, wavelengths1_Keldysh, Efield1_compare*np.sqrt(OpticalIndex_disabled), order, 1., False) #NOTE: Refractive index is accounted in field, not in adiabadicity coefficient
    
    # Dirty: choosing energy or wavelengths. Sorting the results with this choice. 
    if(EnergyAxis):
        wavelengths1_Keldysh_Show = np.divide(2.*np.pi*np.float128(c)*np.float128(hbar)/np.float128(e), wavelengths1_Keldysh) #eV
        wavelengths1_Keldysh_frozen_Show = np.divide(2.*np.pi*c*hbar/e, wavelengths1_Keldysh) #eV
    else:
        wavelengths1_Keldysh_Show = np.multiply(1E9, wavelengths1_Keldysh)
        wavelengths1_Keldysh_frozen_Show = np.multiply(1E9, wavelengths1_Keldysh)
    
    # resorting the results to plot as function of energy if necessary
    lists = sorted(zip(*[wavelengths1_Keldysh_Show, PeakWpi_Gruzdev_continuous]))
    wavelengths1_Keldysh_Show, PeakWpi_Gruzdev_Show = list(zip(*lists)) 
    
    lists_frozen = sorted(zip(*[wavelengths1_Keldysh_frozen_Show, PeakWpi_Gruzdev_frozen]))
    wavelengths1_Keldysh_frozen_Show, PeakWpi_Gruzdev_frozen_Show = list(zip(*lists_frozen)) 
    
    logger.debug(str(np.shape(PeakWpi_Gruzdev_Show))+", "+str(np.shape(wavelengths1_Keldysh_Show)))
    wavelengths_peaks, wpi_peaks, wavelengths_valleys, wpi_valleys=FindMinAndMax(np.array(PeakWpi_Gruzdev_Show), np.array(wavelengths1_Keldysh_Show))
    
    logger.info("W_PI peaks are at: "+str(wavelengths_peaks))
    logger.info("W_PI hills are at: "+str(wavelengths_valleys))
    NumberOfPeaks   = np.shape(wavelengths_peaks)
    NumberOfValleys = np.shape(wavelengths_valleys)
    logger.info("Adding extrema to a list for plotting...")
    if(PlotStarkShift):
        for element in wavelengths_peaks:
            Keldysh_Peaks.append(  [Efield1_compare, element]) #one line per datapoint
        for element in wavelengths_valleys:
            Keldysh_Valleys.append([Efield1_compare, element]) #one line per datapoint
    
    
    plot1 += [plot_field]
  #}}} #end for plotting Fig1a inset.
  
  Keldysh_Peaks=np.array(Keldysh_Peaks).T; Keldysh_Valleys=np.array(Keldysh_Valleys).T
  #print(np.shape(Keldysh_Peaks), np.shape(Keldysh_Valleys))
  
  if(PlotStarkShift):
      if(ShowInset):
        if(PlotValleysInset):
            ax11.scatter(np.multiply(1E-9,Keldysh_Valleys[0,:]), Keldysh_Valleys[1,:], s=0.3, c="blue", marker='*', label="Valleys Keldysh")
        else:
            ax11.scatter(np.multiply(1E-9,Keldysh_Peaks[0,:]), Keldysh_Peaks[1,:], s=0.3, c="blue", marker='*', label="Peaks Keldysh")
  
  ## Here we would like to try adding blue curves from Keldysh (inset) to the main plot. 
      if(PlotValleysInset):
          EnergiesOfValleys=Keldysh_Valleys[1,:]
          FieldsOfValleys=Keldysh_Valleys[0,:]
          # We need to recover the values of Keldysh in continuum of fields for certain wavelengths [the one from EnergiesOfValleys]
 """     
  elif(not PlotValleysInset): #Valleys can be hills now. 
      EnergiesOfValleys=Keldysh_Peaks[1,:]
      FieldsOfValleys=Keldysh_Peaks[0,:]
      # We need to recover the values of Keldysh in continuum of fields for certain wavelengths [the one from EnergiesOfValleys]
      
      #WpiValleys=libKeldysh.GenerateKeldyshDatabase(Egap_SI, meff_Waroquiers_LDA, h*c/(EnergiesOfValleys*e), FieldsOfValleys, order, 1., False)
      
      #print("Number of fields: ", np.shape(FieldsOfValleys), "Number of Wpi:", np.shape(WpiValleys))
      
      #print("Valley min, max: ", np.min(WpiValleys), np.max(WpiValleys))
      
      #ax1.scatter(EnergiesOfValleys, WpiValleys, s=0.3, c="black", marker='*', label="")
  if(Use128bitsFloats):
      ContinuousKeldysh_Wavelengths=np.float128(h)*np.float128(c)/(np.float128(EnergiesOfValleys)*np.float128(e))
  else:
      ContinuousKeldysh_Wavelengths=(h*c)/(EnergiesOfValleys*e)
  WpiValleys=libKeldysh.GenerateKeldyshDatabase(Egap_SI, meff_Waroquiers_LDA, ContinuousKeldysh_Wavelengths, FieldsOfValleys, order, 1., False)
  logger.info("Number of fields: "+str(np.shape(FieldsOfValleys))+"Number of Wpi:"+str(np.shape(WpiValleys)))
  logger.info("Valley min, max: "+str(np.min(WpiValleys))+", "+str(np.max(WpiValleys)))
  ax1.scatter(EnergiesOfValleys, WpiValleys, s=2, c="grey", marker='.', label="")
  ax2.scatter(EnergiesOfValleys, WpiValleys, s=2, c="grey", marker='.', label="")
  #ax1.plot(EnergiesOfValleys, WpiValleys, '.', linewidth=0.1*LineWidthCoeff_rel*LineWidthCoeff_abs, color="grey", label="")
  #ax2.plot(EnergiesOfValleys, WpiValleys, '.', linewidth=0.1*LineWidthCoeff_rel*LineWidthCoeff_abs, color="grey", label="")
  # WARNING: where are the jumps in energy coming from? These energies were captured as minimum of Keldysh function. Re-introducing them into Keldysh is a bit stupid maybe. 
  """
      elif(not PlotValleysInset): #Valleys can be hills now.
          EnergiesOfValleys=Keldysh_Peaks[1,:]
          FieldsOfValleys=Keldysh_Peaks[0,:]
          # We need to recover the values of Keldysh in continuum of fields for certain wavelengths [the one from EnergiesOfValleys]

          #WpiValleys=libKeldysh.GenerateKeldyshDatabase(Egap_SI, meff_Waroquiers_LDA, h*c/(EnergiesOfValleys*e), FieldsOfValleys, order, 1., False)

          #print("Number of fields: ", np.shape(FieldsOfValleys), "Number of Wpi:", np.shape(WpiValleys))

          #print("Valley min, max: ", np.min(WpiValleys), np.max(WpiValleys))

          #ax1.scatter(EnergiesOfValleys, WpiValleys, s=0.3, c="black", marker='*', label="")
      if(Use128bitsFloats):
          ContinuousKeldysh_Wavelengths=np.float128(h)*np.float128(c)/(np.float128(EnergiesOfValleys)*np.float128(e))
      else:
          ContinuousKeldysh_Wavelengths=(h*c)/(EnergiesOfValleys*e)
      WpiValleys=libKeldysh.GenerateKeldyshDatabase(Egap_SI, meff_Waroquiers_LDA, ContinuousKeldysh_Wavelengths, FieldsOfValleys, order, 1., False)
      logger.info("Number of fields: "+str(np.shape(FieldsOfValleys))+"Number of Wpi:"+str(np.shape(WpiValleys)))
      logger.info("Valley min, max: "+str(np.min(WpiValleys))+", "+str(np.max(WpiValleys)))
      ax1.scatter(EnergiesOfValleys, WpiValleys, s=2, c="grey", marker='.', label="")
      ax2.scatter(EnergiesOfValleys, WpiValleys, s=2, c="grey", marker='.', label="")

      #ax1.plot(EnergiesOfValleys, WpiValleys, '.', linewidth=0.1*LineWidthCoeff_rel*LineWidthCoeff_abs, color="grey", label="")
      #ax2.plot(EnergiesOfValleys, WpiValleys, '.', linewidth=0.1*LineWidthCoeff_rel*LineWidthCoeff_abs, color="grey", label="")
      # WARNING: where are the jumps in energy coming from? These energies were captured as minimum of Keldysh function. Re-introducing them into Keldysh is a bit stupid maybe.

  #if(EnergyAxis):
    #ax1.set_xscale('log')
  ax1.set_yscale('log')
  ax1.set_ylabel(r'$w_{PI}$ (m$^{-3}$ s$^{-1}$)')
  
  # Setting the major ticks first (see https://stackoverflow.com/questions/44078409/matplotlib-semi-log-plot-minor-tick-marks-are-gone-when-range-is-large)
  #locmaj = mp.ticker.LogLocator(base=10,numticks=100)
  #ax1.xaxis.set_major_locator(locmaj)
  locmaj = mp.ticker.LogLocator(base=10,numticks=100)
  ax1.yaxis.set_major_locator(locmaj)
  
  # Setting the minor ticks then
  # TODO: solution for minor ticks is here: https://stackoverflow.com/questions/44078409/matplotlib-semi-log-plot-minor-tick-marks-are-gone-when-range-is-large
  #subdiv=np.arange(0.,10.,0.5)
  #locmin = mp.ticker.LogLocator(base=10.0,subs=subdiv,numticks=1000)
  #ax1.xaxis.set_minor_locator(locmin)
  #ax1.xaxis.set_minor_formatter(mp.ticker.NullFormatter())
  subdiv=np.arange(0.,10.,1)
  locmin = mp.ticker.LogLocator(base=10.0,subs=subdiv,numticks=1000)
  ax1.yaxis.set_minor_locator(locmin)
  ax1.yaxis.set_minor_formatter(mp.ticker.NullFormatter())
  
  #ax1.yaxis.set_minor_locator(mp.ticker.LogLocator(base=10.0, subs='auto', numticks=None))
  #ax1.yaxis.set_minor_locator(mp.ticker.LogLocator())
  if(ShowInset):
    ax11.set_xlabel(r'$E_{\mathrm{field}}$ (V/nm)')
    ax11.set_ylabel(r'$\hbar \omega$ (eV)')
  #if(PlotValleysInset):
      #ax11.set_title(r"Min of $\tilde{w}_{PI}$")
  #else:
      #ax11.set_title(r"Max of $\tilde{w}_{PI}$")
  ##ax11.set_xscale('log')
  #ax11.set_yscale('log')
  ax2.set_yscale('symlog', linthresh=1E20)
  ax2.set_ylabel(r'$w_{PI}$ (m$^{-3}$ s$^{-1}$)')
  
  # Setting the major ticks first (see https://stackoverflow.com/questions/44078409/matplotlib-semi-log-plot-minor-tick-marks-are-gone-when-range-is-large)
  #locmaj = mp.ticker.LogLocator(base=10,numticks=100)
  #ax2.xaxis.set_major_locator(locmaj)
  locmaj = mp.ticker.LogLocator(base=10,numticks=100)
  ax2.yaxis.set_major_locator(locmaj)
  
  # Setting the minor ticks then
  # TODO: solution for minor ticks is here: https://stackoverflow.com/questions/44078409/matplotlib-semi-log-plot-minor-tick-marks-are-gone-when-range-is-large
  #subdiv=np.arange(0.,10.,0.5)
  #locmin = mp.ticker.LogLocator(base=10.0,subs=subdiv,numticks=1000)
  #ax2.xaxis.set_minor_locator(locmin)
  #ax2.xaxis.set_minor_formatter(mp.ticker.NullFormatter())
  subdiv=np.arange(0.,10.,1)
  locmin = mp.ticker.LogLocator(base=10.0,subs=subdiv,numticks=1000)
  ax2.yaxis.set_minor_locator(locmin)
  ax2.yaxis.set_minor_formatter(mp.ticker.NullFormatter())
  
  ax2.set_ylabel(r'$w_{PI}$ (m$^{-3}$ s$^{-1}$)')

  # Setting the major ticks first (see https://stackoverflow.com/questions/44078409/matplotlib-semi-log-plot-minor-tick-marks-are-gone-when-range-is-large)
  #locmaj = mp.ticker.LogLocator(base=10,numticks=100)
  #ax2.xaxis.set_major_locator(locmaj)
  locmaj = mp.ticker.LogLocator(base=10,numticks=100)
  ax2.yaxis.set_major_locator(locmaj)

  # Setting the minor ticks then
  # TODO: solution for minor ticks is here: https://stackoverflow.com/questions/44078409/matplotlib-semi-log-plot-minor-tick-marks-are-gone-when-range-is-large
  #subdiv=np.arange(0.,10.,0.5)
  #locmin = mp.ticker.LogLocator(base=10.0,subs=subdiv,numticks=1000)
  #ax2.xaxis.set_minor_locator(locmin)
  #ax2.xaxis.set_minor_formatter(mp.ticker.NullFormatter())
  subdiv=np.arange(0.,10.,1)
  locmin = mp.ticker.LogLocator(base=10.0,subs=subdiv,numticks=1000)
  ax2.yaxis.set_minor_locator(locmin)
  ax2.yaxis.set_minor_formatter(mp.ticker.NullFormatter())

  ax22.set_yscale('symlog', linthresh=1E20)
  
  if(ShowDensities): #{{{
    ax12 = ax1.twinx()
    #plt.title('Final density as function of laser wavelength')
    logger.info(Header+"Total number of datasets: "+str(len(set(Efield1_uniq))))
    
    for Efield1_value in Efield1_uniq: #for each value of Efield1, we make a new search, and plot one new line in the chart
        #PlottingSinglePulses_Efield = libDatabase.FilterDatabase(PlottingSinglePulses, Efield1_value, IndexTable['Efield1']) #select the local Efield
        # TODO: how to merge data with values close to each others? Maybe introduce tolerance of 10% ?




        PlottingSinglePulses_Efield = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, Efield1_value*(1.+FieldTolerance), libPlotting.IndexTable['Efield1'])

        PlottingSinglePulses_Efield = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses_Efield, Efield1_value*(1.-FieldTolerance), libPlotting.IndexTable['Efield1'])
        
        #print PlottingSinglePulses_Efield
        ## prints the current state of selection and stores the values of sweeped parameters
        runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingSinglePulses_Efield)

        
        if(EnergyAxis): #very dirty!
          wavelengths1=np.divide(2.*np.pi*c*hbar/e, wavelengths1)
        
        lists = sorted(zip(*[wavelengths1, Ne_final, ElectronEnergy_max])) #NOTE: Ne_max makes more sense than Ne_final
        wavelengths1_sorted, Ne_final_sorted, Energy_max_sorted = list(zip(*lists))
        if(not EnergyAxis): 
            wavelengths1_sorted = np.multiply(1e9, wavelengths1_sorted)
        
        ax12.set_ylabel(r'Final density $n_{\mathrm{exc}}$ (m$^{-3}$)')
        #print Efield1_value
        # TODO: PUT here a select case which attributes a new color to a new field interval
        # OR USE A DICTIONNARY. 
        
        # SAME CODE, just more clever. 
        for Efield1_compare in FieldOfInterest:
            if(abs(Efield1_value - Efield1_compare)/Efield1_value < FieldTolerance):
                #if(not EnergyAxis):
                plot_field, = ax12.plot(wavelengths1_sorted, Ne_final_sorted, '^', color=Colors[Efield1_compare], label=r'$E_{peak}=$'+str(round(Efield1_compare*1E-9))+' V/nm')

                #plot_field, = ax12.plot(wavelengths1_sorted, Ne_final_sorted, '^', label=r'$E_{peak}=$'+str(round(Efield1_compare*1E-9))+' V/nm')
                plot1 += [plot_field]
                #, I='+str(libKeldysh.FieldToIntensity(Efield1_value))+r' W/m$^2$')
                break
        
        del PlottingSinglePulses_Efield
  #}}}  
  labels1 = [l.get_label() for l in plot1]
  #ax1.legend(plot1, labels1, loc='best') #outside right, middle
  #ax1.legend(loc=3, borderaxespad=0., bbox_to_anchor=(1.05, 0.05))
  if(not EnergyAxis): 
    ax1.set_xlabel('Laser wavelength (nm)')
    ax1.set_xlim([wavelengths1_min, wavelengths1_max])
    if(ShowInset): 
        ax11.set_ylim([wavelengths1_min,wavelengths1_max])
    ax2.set_xlabel('Laser wavelength (nm)')
    ax2.set_xlim([wavelengths1_min, wavelengths1_max])
    ax22.set_xlim([wavelengths1_min,wavelengths1_max])
  else: 
    ax1.set_xlabel('Photon energy (eV)')
    ax1.set_xlim([h*c/wavelengths1_max/e, h*c/wavelengths1_min/e])
    if(ShowInset): 
        #ax11.set_ylim([h*c/wavelengths1_max/e, h*c/wavelengths1_min/e])
        ax11.set_ylim([h*c/wavelengths1_max/e, 1.5])
    ax2.set_xlabel('Photon energy (eV)')
    ax2.set_xlim([h*c/wavelengths1_max/e, h*c/wavelengths1_min/e])
    ax22.set_xlim([h*c/wavelengths1_max/e, h*c/wavelengths1_min/e])
  if(AccountForOpticalIndex): 
      ax1.set_ylim([1E35,1E43])
      ax2.set_ylim([1E35,1E43])
      ax22.set_ylim([1E35,1E43])
  else: 
      ax1.set_ylim([1E32,1E43])
      ax2.set_ylim([1E32,1E43])
      ax22.set_ylim([1E32,1E43])
  #ax1.set_xscale('log')
  #ax1.set_yscale('log')
  #ax1.grid()
  
  ResonantWavelength = h*c/Egap_SI
  logger.info("** Adding vertical bar for the resonant transitions at: "+str(round(1e9*ResonantWavelength))+"nm.")
  
  for order in np.arange(1,DrawOrders+1):
    if(not EnergyAxis):
      ax1.axvline(x=1e9*order*ResonantWavelength, color='k', linestyle='--', linewidth=0.5)
      #ax1.axvline(x=1e9*2*ResonantWavelength, color='k', linestyle='--', linewidth=1)
      #ax1.axvline(x=1e9*3*ResonantWavelength, color='k', linestyle='--', linewidth=1)
      #ax1.axvline(x=1e9*4*ResonantWavelength, color='k', linestyle='--', linewidth=1)
      
      ax2.axvline(x=1e9*order*ResonantWavelength, color='k', linestyle='--', linewidth=0.5)
      #ax2.axvline(x=1e9*2*ResonantWavelength, color='k', linestyle='--', linewidth=1)
      #ax2.axvline(x=1e9*3*ResonantWavelength, color='k', linestyle='--', linewidth=1)
      #ax2.axvline(x=1e9*4*ResonantWavelength, color='k', linestyle='--', linewidth=1)
      
      ax22.axvline(x=1e9*order*ResonantWavelength, color='k', linestyle='--', linewidth=1)
      #ax22.axvline(x=1e9*2*ResonantWavelength, color='k', linestyle='--', linewidth=1)
      #ax22.axvline(x=1e9*3*ResonantWavelength, color='k', linestyle='--', linewidth=1)
      #ax22.axvline(x=1e9*4*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    else:
      ax1.axvline(x=Egap_SI/e/float(order), color='k', linestyle='--', linewidth=1)
      #ax1.axvline(x=1./2.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
      #ax1.axvline(x=1./3.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
      #ax1.axvline(x=1./4.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
      if(ShowInset): 
        ax11.axhline(y=Egap_SI/e/float(order), color='k', linestyle='--', linewidth=1)
        #ax11.axhline(y=1./2.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
        #ax11.axhline(y=1./3.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
        #ax11.axhline(y=1./4.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
      
      ax2.axvline(x=Egap_SI/e/float(order), color='k', linestyle='--', linewidth=1)
      #ax2.axvline(x=1./2.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
      #ax2.axvline(x=1./3.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
      #ax2.axvline(x=1./4.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
      
      ax22.axvline(x=Egap_SI/e/float(order), color='k', linestyle='--', linewidth=1)
      #ax22.axvline(x=1./2.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
      #ax22.axvline(x=1./3.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
      #ax22.axvline(x=1./4.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
  
  ## Saving this first double figure
  plt.tight_layout()
  filename = 'Metaplot-SinglePulse-SingleColor-EffectOfWavelength_Tol10-Keldysh-'+str(int(pulseDuration*1E15))
  if(PlotValleysInset):
      filename+="-Valleys"
  else:
      filename+="-Hills"
  if(AccountForOpticalIndex): 
      filename+="-WithOpticalIndex"
  else:
      filename+="-NoOpticalIndex"
  if(PlotSaturation):
      filename+="-WithSaturation"
  plt.savefig(filename+'.eps')
  plt.savefig(filename+'.png')
  
  ## Imports the optical data from Palik database
  # Returns the spectrum (in Re(eps), Im(eps)) of a material on the given grid of wavelengths
  PulseDurationLocal=pulseDuration

  ## Computes the free carrier energy contained in the volume <SampleVolume>.  
  def EnergyFromFreeCarrierAbsorption(wavelengths1_sorted, Efield, Ne_final_sorted, nu=6e-15**-1):
                                      #Egap_SI = 2.56*e, meff = 0.2226, SampleVolume=2 * 5E28**-1): #{{{
    # Since parameters are not welcome in function adjustment, I need to pass some as global variables
    global Egap_SI, meff, SampleVolume, pulseDuration

    # dEnergy/dt = alpha * intensity + ( n har omega - Egap ) x sigma_n I^n
    # where we approach sigma_n I^n by pulse-averaged w_PI 
    # Energy ~ alpha * intensity + ( n har omega - Egap ) x w_PI * tau
    omega   = np.divide(2. * np.pi * c, wavelengths1_sorted)
    #print Header+"omega: "+str(omega) #Valid
    # nu = 6e-15**-1 #1.1 fs: Solowski-Tinten adjustment. Can also use the one from Shunsuke Sato, PRB 2015 (Numerical pump probe)

    logger.info("nu: "+str(nu)) #Valid
    #epsilon_Si0 = 13.64 + 0.048j #TODO: get the optical data in extended spp theory for Si
    epsilon_Si0 = libOpticalData.ImportPalikDatabase_epsilon_fromNK(wavelengths1_sorted) ## TODO: replace by LDA !
    logger.debug("Dielectric permittivity of Si at 800 nm : "+str(libOpticalData.ImportPalikDatabase_epsilon_fromNK(800e-9)))
    omega_p2    = np.multiply(Ne_final_sorted, e**2 / m_e / meff / epsilon_0)
    #print "omega_p: "+str(np.sqrt(omega_p2)) # Valid
    # NOTE: THIS IS A CHOICE. WE JUST HAVE TO BE CONSISTENT ON WHICH PROCESS WE WISH FOR. 
    Epsilon = 1E0 - np.divide(omega_p2, np.power(omega, 2)) * np.divide(1.,1. + np.divide(1.j * nu, omega))
    #Epsilon = epsilon_Si0 - np.divide(omega_p2, np.power(omega, 2)) * np.divide(1.,1. + np.divide(1.j * nu, omega))
    #NumberOfPhotons = Egap_SI / (hbar * omega)
    NumberOfPhotons = np.divide(Egap_SI/hbar, omega)
    logger.info(NumberOfPhotons)
    NumberOfPhotons_round = np.ceil(NumberOfPhotons)
    logger.debug(Header+"sqrt(Epsilon): "+str(np.min(np.sqrt(Epsilon)))+", "+str(np.max(np.sqrt(Epsilon)))) #Valid
    alpha                = np.multiply(2.* omega/c, np.imag(np.sqrt(Epsilon))) #valid
    EnergyExcess         = np.multiply(NumberOfPhotons_round * hbar * omega, Ne_final_sorted) * SampleVolume / e  #result in eV
    EnergyBremsstrahlung = np.multiply( omega * epsilon_0 * (Epsilon.imag), np.multiply(Efield, Efield)) * SampleVolume / e  * pulseDuration
    EnergyBremsstrahlung_approx =  EnergyBremsstrahlung #NOTE: not used for now. Another collision freq? 
    #EnergyBremsstrahlung_approx = np.multiply( 4.*np.pi / c * np.divide(Epsilon.imag, (np.sqrt(Epsilon)).real), )) ) * SampleVolume / e
    logger.info("FcrEnergy (eV): "+str(np.min(EnergyBremsstrahlung))+", "+str(np.max(EnergyBremsstrahlung)))
    logger.info("Excess energy (n hbar omega) (eV): "+str(np.min(EnergyExcess))+", "+str(np.max(EnergyExcess)))
    EnergyFreeCarrier    = EnergyExcess + EnergyBremsstrahlung
    EnergyFreeCarrier_approx = EnergyExcess + EnergyBremsstrahlung_approx
    return EnergyFreeCarrier, EnergyFreeCarrier_approx
  #}}}
  
  
  
  
  ## TODO: TDDFT vs Drude should become another figure. We don't need to have them side by side. 
  plt.figure()
  ax3 = plt.subplot(111)
  #plt.title('Irreversibility factor')
  #plt.plot(libKeldysh.FieldToIntensity(Efield1), irreversibility, 'o', label='Irreversibility')
  for Efield1_value in Efield1_uniq: #for each value of Efield1, we make a new search, and plot one new line in the chart
    #PlottingSinglePulses_Efield = libDatabase.FilterDatabase(PlottingSinglePulses, Efield1_value, IndexTable['Efield1']) #select the local Efield




    PlottingSinglePulses_Efield = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, Efield1_value*(1.+FieldTolerance), libPlotting.IndexTable['Efield1'])

    PlottingSinglePulses_Efield = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses_Efield, Efield1_value*(1.-FieldTolerance), libPlotting.IndexTable['Efield1'])
    #print Header+"Dataset for the field: "+str(Efield1_value*1e-9)+' V/nm'
    #print PlottingSinglePulses_Efield
    # TODO: as fields are all slightly different, we should project on the neareast neighboring value, up to 1E8 V/m.
    
    ## prints the current state of selection and stores the values of sweeped parameter

    runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingSinglePulses_Efield)

    
    # 
    if(EnergyAxis):
        wavelengths1_Show = np.divide(2.*np.pi*c*hbar/e, wavelengths1) #eV
    else: 
        wavelengths1_Show = np.multiply(1e9, wavelengths1)
    
    lists = sorted(zip(*[wavelengths1_Show, wavelengths1, Ne_final, W_pi, ElectronEnergy_max/NumberOfAtoms])) #NOTE: Ne_max makes more sense than Ne_final
    wavelengths1_sorted_show, wavelengths1_sorted, Ne_final_sorted, W_pi_sorted, Energy_max_sorted = list(zip(*lists))
    #print Efield1_value
    
    for Efield1_compare in FieldOfInterest: 
        if(abs(Efield1_value - Efield1_compare)/Efield1_value < FieldTolerance):
            
            ax3.plot(wavelengths1_sorted_show, Energy_max_sorted, '-', linewidth=LineWidthCoeff_abs*LineWidths[Efield1_compare], color=Colors[Efield1_compare])
            #ax3.plot(wavelengths1_sorted_show, Energy_max_sorted, '-')
            
    del PlottingSinglePulses_Efield
    
  #plt.xlabel(r'Intensity (W/m$^2$)')
  if(not EnergyAxis):
    ax3.set_xlabel('Laser wavelength (nm)')
  else: 
    ax3.set_xlabel('Photon energy (eV)')
  #plt.ylabel('Irreversibility factor')
  ax3.set_ylabel('Energy absorbed by electrons (eV)')
  #plt.ylim((0,1))
  if(not EnergyAxis): 
    ax3.set_xlim([wavelengths1_min, wavelengths1_max])
    plt.xlim((wavelengths1_min, wavelengths1_max))#ax2.set_xscale('log')
  else: 
    ax3.set_xlim([h*c/wavelengths1_max/e, h*c/wavelengths1_min/e])
    plt.xlim((h*c/wavelengths1_max/e, h*c/wavelengths1_min/e))#ax2.set_xscale('log')
  
  ## TODO: this can be removed, since irreversibility is gauge-dependent. 
  #if(PlotIrreversibility): #{{{
    #ax22 = ax2.twinx()
    #for Efield1_value in Efield1_uniq: #for each value of Efield1, we make a new search, and plot one new line in the chart
        ##PlottingSinglePulses_Efield = libDatabase.FilterDatabase(PlottingSinglePulses, Efield1_value, libPlotting.IndexTable['Efield1']) #select the local Efield
        #PlottingSinglePulses_Efield = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, Efield1_value*(1.+FieldTolerance), libPlotting.IndexTable['Efield1'])
        #PlottingSinglePulses_Efield = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses_Efield, Efield1_value*(1.-FieldTolerance), libPlotting.IndexTable['Efield1'])
        ##print Header+"Dataset for the field: "+str(Efield1_value*1e-9)+' V/nm'
        ##print PlottingSinglePulses_Efield
        ## TODO: as fields are all slighly different, we should project on the neareast neighboring value, up to 1E8 V/m. 
        
        ### prints the current state of selection and stores the values of sweeped parameters
        #runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, ElectronEnergy_max, Energy_max, W_pi, PulseType, RiseTime, FinalTime , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingSinglePulses_Efield)
        
        #lists = sorted(itertools.izip(*[wavelengths1, Ne_final, Energy_max, irreversibility])) #NOTE: Ne_max makes more sense than Ne_final
        #wavelengths1_sorted, Ne_final_sorted, Energy_max_sorted, irreversibility_sorted = list(itertools.izip(*lists))
        #wavelengths1_sorted_show = np.multiply(1e9, wavelengths1_sorted)
        
        #for Efield_compare in FieldOfInterest:
            #if(abs(Efield1_value - Efield1_compare)/Efield1_compare < FieldTolerance):
                #ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+'+Colors[Efield1_compare], label='Irreversibility')
        
        ###print Efield1_value
        ### TODO: PUT here a select case which attributes a new color to a new field interval
        ###if(1e8*(1.-FieldTolerance) < Efield1_value < 1e8*(1.+FieldTolerance)):
            ###ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+r', label='Irreversibility')
            ###ax22.plot(wavelengths1_sorted_show, libKeldysh.gammaKeldysh(Egap_SI, 0.228, Efield1_value, wavelengths1_sorted), '--r', label='Adiabaticity')
        ##if(1e9*(1.-FieldTolerance) < Efield1_value < 1e9*(1.+FieldTolerance)):
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+b', label='Irreversibility')
            ###ax22.plot(wavelengths1_sorted_show, libKeldysh.gammaKeldysh(Egap_SI, 0.228, Efield1_value, wavelengths1_sorted), '--b', label='Adiabaticity')
        ##if(5e9*(1.-FieldTolerance) < Efield1_value < 5e9*(1.+FieldTolerance)):
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+g', label='Irreversibility')
            ###ax22.plot(wavelengths1_sorted_show, libKeldysh.gammaKeldysh(Egap_SI, 0.228, Efield1_value, wavelengths1_sorted), '--g', label='Adiabaticity')
        ##if(1e10*(1.-FieldTolerance) < Efield1_value < 1e10*(1.+FieldTolerance)):
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+k', label='Irreversibility')
            ###ax22.plot(wavelengths1_sorted_show, libKeldysh.gammaKeldysh(Egap_SI, 0.228, Efield1_value, wavelengths1_sorted), '--k', label='Adiabaticity')
        #del PlottingSinglePulses_Efield
    
    
    #plt.ylabel(r'Irreversibility $\zeta$')
    #ax22.set_yscale('log')
  ##}}}
  

  if(PlotFreeCarrierAbsorption):
      SampleVolume = 2 * 5E28 ** -1 #used in Drude function
      logger.info("Now plotting the free carrier energy absorption...")
      for Efield1_value in Efield1_uniq: #for each value of Efield1, we make a new search, and plot one new line in the chart
        #PlottingSinglePulses_Efield = libDatabase.FilterDatabase(PlottingSinglePulses, Efield1_value, libPlotting.IndexTable['Efield1']) #select the local Efield




        PlottingSinglePulses_Efield = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, Efield1_value*(1.+FieldTolerance), libPlotting.IndexTable['Efield1'])

        PlottingSinglePulses_Efield = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses_Efield, Efield1_value*(1.-FieldTolerance), libPlotting.IndexTable['Efield1'])
        #print Header+"Dataset for the field: "+str(Efield1_value*1e-9)+' V/nm'
        #print PlottingSinglePulses_Efield
        # TODO: as fields are all slighly different, we should project on the neareast neighboring value, up to 1E8 V/m. 
        
        ## prints the current state of selection and stores the values of sweeped parameters



        runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, \
        pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, \
        Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, \
        TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, \
        HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, \
        HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, \
        W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingSinglePulses_Efield)

        
        if(EnergyAxis):
            wavelengths1_Show = np.divide(2.*np.pi*c*hbar/e, wavelengths1) #eV
        else: 
            wavelengths1_Show = np.multiply(1e9, wavelengths1)
    
        lists = sorted(zip(*[wavelengths1_Show, wavelengths1, Ne_final, W_pi, ElectronEnergy_max])) #NOTE: Ne_max makes more sense than Ne_final
        wavelengths1_sorted_show, wavelengths1_sorted, Ne_final_sorted, W_pi_sorted, Energy_max_sorted = list(zip(*lists))
        #print Efield1_value
    
        
        #lists = sorted(itertools.izip(*[wavelengths1, Ne_final, W_pi, Energy_max, irreversibility])) #NOTE: Ne_max makes more sense than Ne_final
        #wavelengths1_sorted, Ne_final_sorted, W_pi_sorted, Energy_max_sorted, irreversibility_sorted = list(itertools.izip(*lists))
        #wavelengths1_sorted_show = np.multiply(1e9, wavelengths1_sorted)
        nu_t=[] #1 value per FieldOfInterest
        idx=0
        for Efield1_compare in FieldOfInterest:
            if(abs(Efield1_value - Efield1_compare) / Efield1_value < FieldTolerance):
                """
                EnergyAnalytical, EnergyAnalytical_approx = EnergyFromFreeCarrierAbsorption(pulses1, wavelengths1_sorted, Efield1_value, Ne_final_sorted, W_pi_sorted, Egap_SI, EffectiveMass[Efield1_compare])
                                
                ax3.plot(wavelengths1_sorted_show, EnergyAnalytical/NumberOfAtoms, '--', color=Colors[Efield1_compare], label='Free carrier absorption')
                #ax2.plot(wavelengths1_sorted_show, EnergyAnalytical, '--', label='Free carrier absorption')
            
        ##if(1e8*(1.-FieldTolerance) < Efield1_value < 1e8*(1.+FieldTolerance)):
            ##print Header+"Found field 0.1 V/nm"
            ##EnergyAnalytical, EnergyAnalytical_approx = EnergyFromFreeCarrierAbsorption(pulses1, wavelengths1_sorted, Efield1_value, Ne_final_sorted, W_pi_sorted, Egap_SI)
            ##ax2.plot(wavelengths1_sorted_show, EnergyAnalytical, '--+r', label='Free carrier absorption')
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+r', label='Irreversibility')
            ##ax22.plot(wavelengths1_sorted_show, libKeldysh.gammaKeldysh(Egap_SI, 0.228, Efield1_value, wavelengths1_sorted), '--r', label='Adiabaticity')
        #if(1e9*(1.-FieldTolerance) < Efield1_value < 1e9*(1.+FieldTolerance)):
            #print Header+"Found field 1.0 V/nm"
            #EnergyAnalytical, EnergyAnalytical_approx = EnergyFromFreeCarrierAbsorption(pulses1, wavelengths1_sorted, Efield1_value, Ne_final_sorted, W_pi_sorted, Egap_SI)
            #ax2.plot(wavelengths1_sorted_show, EnergyAnalytical, '--+b', label='Free carrier absorption')
            ##ax2.plot(wavelengths1_sorted_show, EnergyAnalytical_approx, '--b', label='Free carrier absorption')
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+b', label='Irreversibility')
            ##ax22.plot(wavelengths1_sorted_show, libKeldysh.gammaKeldysh(Egap_SI, 0.228, Efield1_value, wavelengths1_sorted), '--b', label='Adiabaticity')
        #if(5e9*(1.-FieldTolerance) < Efield1_value < 5e9*(1.+FieldTolerance)):
            #print Header+"Found field 5 V/nm"
            #EnergyAnalytical, EnergyAnalytical_approx = EnergyFromFreeCarrierAbsorption(pulses1, wavelengths1_sorted, Efield1_value, Ne_final_sorted, W_pi_sorted, Egap_SI)
            #ax2.plot(wavelengths1_sorted_show, EnergyAnalytical, '--+g', label='Free carrier absorption')
            ##ax2.plot(wavelengths1_sorted_show, EnergyAnalytical_approx, '--g', label='Free carrier absorption')
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+g', label='Irreversibility')
            ##ax22.plot(wavelengths1_sorted_show, libKeldysh.gammaKeldysh(Egap_SI, 0.228, Efield1_value, wavelengths1_sorted), '--g', label='Adiabaticity')
        #if(1e10*(1.-FieldTolerance) < Efield1_value < 1e10*(1.+FieldTolerance)):
            #print Header+"Found field 10 V/nm"
            #EnergyAnalytical, EnergyAnalytical_approx = EnergyFromFreeCarrierAbsorption(pulses1, wavelengths1_sorted, Efield1_value, Ne_final_sorted, W_pi_sorted, Egap_SI)
            #ax2.plot(wavelengths1_sorted_show, EnergyAnalytical, '--+k', label='Free carrier absorption')
            ##ax2.plot(wavelengths1_sorted_show, EnergyAnalytical_approx, '--k', label='Free carrier absorption')
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+k', label='Irreversibility')
            ##ax22.plot(wavelengths1_sorted_show, libKeldysh.gammaKeldysh(Egap_SI, 0.228, Efield1_value, wavelengths1_sorted), '--k', label='Adiabaticity')
                """

                # Let's try to make a self-consistent collision frequency
                #try:
                if(True):
                    logger.warning("Fitting with the following parameters: ")
                    logger.warning(str(pulses1[idx]))
                    logger.warning(str(Efield1_compare))
                    logger.warning(str(Ne_final_sorted[idx]))
                    logger.warning(str(W_pi_sorted[idx]))
                    logger.warning(str(Egap_SI))
                    logger.warning(str(EffectiveMass[Efield1_compare]))
                    logger.warning(str(SampleVolume))
                    logger.warning(str(10e-15**-1))
                    logger.debug("Wavelengths1_sorted: "+str(len(wavelengths1_sorted))+str(wavelengths1_sorted))
                    logger.debug("Energy_max_sorted: "+str(len(Energy_max_sorted))+str(Energy_max_sorted))
                    p0=[Efield1_compare, Ne_final_sorted[idx], 1E-15**-1]

                    # Determine the number of parameters expected by EnergyFromFreeCarrierAbsorption function
                    num_params = 3

                    # Set initial parameter values
                    lower_bounds = np.multiply(0.999, [Efield1_compare, Ne_final_sorted[idx], 10e-15**-1])
                    upper_bounds = np.multiply(1.001, [Efield1_compare, Ne_final_sorted[idx], 0.1E-15 ** -1])
                    bounds = (lower_bounds[:num_params], upper_bounds[:num_params])

                    FitDrudeModel_optparameters, FitDrudeModel_covariance = curve_fit(EnergyFromFreeCarrierAbsorption,
                                                         xdata=wavelengths1_sorted,
                                                         ydata=Energy_max_sorted,
                                                         p0=[Efield1_compare, Ne_final_sorted[idx], 1E-15**-1],
                                                         bounds=(np.multiply(0.999,[Efield1_compare, Ne_final_sorted, 10e-15**-1]),
                                                                 np.multiply(1.001,[Efield1_compare, Ne_final_sorted, 0.1E-15**-1])))

                    # def EnergyFromFreeCarrierAbsorption(wavelengths1_sorted, Efield, Ne_final_sorted, nu=6e-15**-1):
                    # We have to minimize ElectronEnergy - Drude model
                    # print Header+"FITTING DRUDE MODEL."
                    # print np.divide(1.,FitDrudeModel_parameters)
                    nu_t = nu_t.append([FitDrudeModel_parameters[8]]) #capture last value
                    logger.warning("Fitting by Drude: "+str(nu_t))

                    #for nuLocal in nu_t: #but then it will make 1 fit per plotted field. Maybe even better?
                    EnergyAnalytical, EnergyAnalytical_approx = EnergyFromFreeCarrierAbsorption(wavelengths1_sorted, *FitDrudeModel_optparameters)

                    #except:
                    """else:
                        nu=6E-15**-1 #value by default which fits rather well the results for all fields
                        logger.warning("Couldn't fit the collision frequency for field"+str(Efield1_compare)+". Taking value by default: "+str(nu**-1)+" fs")
                        EnergyAnalytical, EnergyAnalytical_approx = EnergyFromFreeCarrierAbsorption(wavelengths1_sorted,
                                                                                                    pulses1, Efield1_value,
                                                                                                    Ne_final_sorted,
                                                                                                    W_pi_sorted,
                                                                                                    Egap_SI,
                                                                                                    EffectiveMass[Efield1_compare],
                                                                                                    SampleVolume, nu)"""
                    if(PlotFreeCarrierAbsorption):
                        ax3.plot(wavelengths1_sorted_show, EnergyAnalytical/NumberOfAtoms, '--', color=Colors[Efield1_compare], label='Free carrier absorption')
                        idx += 1 #This is an index to follow desired field curve
                    #ax2.plot(wavelengths1_sorted_show, EnergyAnalytical, '--', label='Free carrier absorption')
# Modules
        del PlottingSinglePulses_Efield
    
    #plt.ylabel(r'Irreversibility $\zeta$')
    #ax22.set_yscale('log')
      #ax22 = ax2.twinx()
  #plt.ylabel(r'Adiabaticity $\gamma$, Irreversibility $\zeta$')
  #ax22.set_ylim((0,1))
  ax3.set_xscale("log")
  ax3.set_yscale("symlog", linthresh=1e-9)
  ax3.yaxis.set_minor_locator(mp.ticker.LogLocator(base=10.0, subs='auto', numticks=None))
  #ax2.set_ylim((1E-9,5E2))
  #ax22.grid()
  if(not EnergyAxis): 
    plt.xlim((200.,2100.))
  else: 
    plt.xlim((0.3,10.))
  
  
  #ResonantWavelength = h*c/Egap_SI
  logger.info("Adding vertical bar for the resonant transitions at: "+str(round(1e9*ResonantWavelength))+"nm.")
  
  for order in np.arange(1,DrawOrders+1):
    if(not EnergyAxis):
      ax3.axvline(x=1e9*order*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    #ax3.axvline(x=1e9*2*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    #ax3.axvline(x=1e9*3*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    #ax3.axvline(x=1e9*4*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    else:
      ax3.axvline(x=Egap_SI/e/float(order), color='k', linestyle='--', linewidth=1)
    #ax3.axvline(x=1./2.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    #ax3.axvline(x=1./3.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    #ax3.axvline(x=1./4.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
  #ax2.legend(loc=2, borderaxespad=0.)
  #ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
  #ax22.legend(loc=1, borderaxespad=0.)

  # Saving the plot
  plt.tight_layout
  filename = 'Metaplot-SinglePulse-SingleColor-EffectOfWavelength_Tol10-Drude-'+str(int(pulseDuration*1E15))
  plt.savefig(filename+'.eps')
  plt.savefig(filename+'.png')
  if(ShowFinalPlot==1):
    plt.show()

  del ax1, PlottingSinglePulses
  #}}}

## Plots (Wpi | Energy) (wavelength) (1D plot) for the special case of 2-color irradiations. Fields are selected by introducing <FieldToTolerance> percents tolerance on the electric field amplitude. 
# This enables to analyze the role of laser wavelength on the final excited electron density, excitation rate and laser-induced energy increase. Provides predictions on the Stark effect. 
# Energy is compared with analytical Drude model. 
def Bicolor_plotEffectOfWavelength_Tolerance(PlottingTable, wavelength1, FieldTolerance=0.05, Field1=1E0, Field2max=1E9, pulseDuration=20e-15, CEP1=0., PulseType=3, RiseTimeRequested=2, ShowFinalPlot=0, Egap_SI=2.56*e, Egap_SI_indirect=0.51*e): #{{{
  ShowDensities            =False
  PlotIrreversibility      =False
  PlotFreeCarrierAbsorption=False
  EnergyAxis               =False #True: Xaxis is in eV. False: wavelength in nm. 
  
  CutOffWavelength=2.1e-6
  
  CEP2_request=-np.pi/2. #np.pi/2.
  
  meff_Waroquiers_LDA = 0.2226 #effective mass of electron [Waroquiers, D. et al, Physical Review B, 2013, 87, 075121]
  #meff = 0.14 #0.14: fits for Si in MPI regime at 800 nm.
  #Field dependent fitting of effective mass
  meff_LowField = meff_Waroquiers_LDA
  meff_HighField = meff_Waroquiers_LDA #0.14: fits for Si in MPI regime at 800 nm.
  meff_VeryHighField = meff_Waroquiers_LDA
  
  meff_Tunnel = meff_Waroquiers_LDA #fits well if using INDIRECT band gap for 800 nm
  print("#====== Plot: DOUBLE PULSE, EFFECT OF THE WAVELENGTH, SAME CEP, SAME POLARIZATIONS: final density as function of field ======")
  #Filter: double pulses, same wavelengths, same polarization
  
  LineWidthCoeff_abs = 0.5 #controls absolute thickness of all curves
  LineWidthCoeff_rel = 1.0 #controls thickness of Keldysh-Gruzdev curves relatively to TDDFT
  antialias=True
  
  #FieldOfInterest = [5e8, 1e9, 2e9, 3E9, 4E9, 5e9, 1e10]
  FieldOfInterest = [2.5E8, 5E8, 1E9] #, 5E9] #[0E0, 1.0E9] #, 1e10] # NOTE: this should be used only for Efield2
  NumberOfColors=len(FieldOfInterest)
  ColorList=[plt.cm.seismic(i) for i in np.linspace(0, 1, NumberOfColors)] 
  # Not bad! RdYlBu, coolwarm, seismic
  # Too soft: Copper, Spectral
  # Not perfect. RdBu. RdYlGn.
  # Too dynamic: bwr
  # To test: , brg, rainbow
  #Colors={5e8: 'c', 7.5e8: 'r', 1e9: 'b', 1.5e9: 'g', 2e9: 'y', 3e9: 'm'} #, 1e10: 'k'}
  Colors={FieldOfInterest[i]:ColorList[i] for i in range(0,NumberOfColors)}
  LineWidths={FieldOfInterest[i]:np.arange(1.5,4,NumberOfColors) for i in range(0,NumberOfColors)}
  #EffectiveMass={1e9: meff_LowField, 2e9: meff_LowField, 5e9: meff_HighField, 1e10: meff_VeryHighField}
  EffectiveMass={FieldOfInterest[i]:meff_HighField for i in range(0, NumberOfColors)}
  # Simple pulse, single color
  
  # TODO: Filtering function
  PlottingDoublePulses = PlottingTable #get the full database
  
  # The wavelength1 parameter enables to choose which wavelength1 we want to observe the Stark effect for. Of course, we shall not discernate laser1 and laser2 in the data. But how? Lets discernate first. And invert later. 
  PlottingDoublePulses = libDatabase.FilterDatabase(PlottingDoublePulses, wavelength1, libPlotting.IndexTable['wavelengths1'])
  
  PlottingDoublePulses = libDatabase.FilterDatabaseGreaterThan(PlottingDoublePulses, 1E0, libPlotting.IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. # NOTE: works
  
  # NOTE: 0 should be kept here, since we use non-strict inequality
  PlottingDoublePulses = libDatabase.FilterDatabaseGreaterThan(PlottingDoublePulses, Field1*(1-FieldTolerance), libPlotting.IndexTable['Efield1']) #select fields that are close to selected Efield1
  PlottingDoublePulses = libDatabase.FilterDatabaseLowerThan(PlottingDoublePulses, Field1*(1+FieldTolerance), libPlotting.IndexTable['Efield1']) #select fields that are close to selected Efield1
  
  # We keep Efield2==0 possible. 
  PlottingDoublePulses = libDatabase.FilterDatabaseLowerThan(PlottingDoublePulses, Field2max, libPlotting.IndexTable['Efield2']) #we want to include cases where Efield2==0, and Efield2 >> 1. However, we want to put a maximum field for Efield2.

  
  PlottingDoublePulses = libDatabase.FilterDatabase(PlottingDoublePulses, 0., libPlotting.IndexTable['CrossPolarization']) #removes fields where CrossPolarization is different from 0. # TODO: CrossPolarization=1 case would be interesting for Stark effect.
  
  PlottingDoublePulses = libDatabase.FilterDatabase(PlottingDoublePulses, PulseType, libPlotting.IndexTable['PulseType']) #keeps only pulses of type StH (3)
  
  PlottingDoublePulses = libDatabase.FilterDatabase(PlottingDoublePulses, pulseDuration, libPlotting.IndexTable['pulses1']) #keeps only pulses of selected duration in argument
  
  if(RiseTimeRequested >= 1): #select by RiseTime only if it is defined
    PlottingDoublePulses = libDatabase.FilterDatabase(PlottingDoublePulses, RiseTimeRequested, libPlotting.IndexTable['RiseTime'])
  
  PlottingDoublePulses = libDatabase.FilterDatabase(PlottingDoublePulses, CEP1, libPlotting.IndexTable['CEP1']) #keeps only CEP = 0.
  
  if(CEP2_request >= 0):
      PlottingDoublePulses = libDatabase.FilterDatabaseLowerThan(PlottingDoublePulses,   (CEP2_request)*(1.1), libPlotting.IndexTable['CEP2'])
      PlottingDoublePulses = libDatabase.FilterDatabaseGreaterThan(PlottingDoublePulses, (CEP2_request)*(0.9), libPlotting.IndexTable['CEP2']) #keeps only CEP = 0.
  else:
      PlottingDoublePulses = libDatabase.FilterDatabaseGreaterThan(PlottingDoublePulses, (CEP2_request)*(1.1), libPlotting.IndexTable['CEP2'])


      PlottingDoublePulses = libDatabase.FilterDatabaseLowerThan(PlottingDoublePulses,   (CEP2_request)*(0.9), libPlotting.IndexTable['CEP2']) #keeps only CEP = 0.
  # TODO: we shall clean CEP2 as well, as it must influence the results. 
  
  PlottingDoublePulses = libDatabase.FilterDatabaseLowerThan(PlottingDoublePulses, CutOffWavelength, libPlotting.IndexTable['wavelengths2']) #introducing a cutoff at 2.1 um in selection.

  ## TODO: check if CEP could count at smaller wavelength. 

  logger.info("Available data: ")
  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = libPlotting.CountTheNumberOfLaserParameters(PlottingDoublePulses)

  #if(len(CEP1_uniq) == 1): #{{{ #single color simulations
  # Complexity level: 2
  
  #print("Number of Efield2 uniq: ")
  #print np.shape(Efield2_uniq)
  
  logger.info(CEP2_uniq)
  #exit()
  
  plt.figure(figsize=(12,6)) #new figure
  plotA=121; plotB=122
  
  ax1 = plt.subplot(111)
  
  #ax1.set_prop_cycle('color',Colors)
    
  plot1 = [] #np.empty((0,0))
  for Efield2_compare in FieldOfInterest: #{{{
    PlottingDoublePulses_Efield = PlottingDoublePulses
    PlottingDoublePulses_Efield = libDatabase.FilterDatabaseLowerThan(PlottingDoublePulses_Efield, Efield2_compare*(1.+FieldTolerance), libPlotting.IndexTable['Efield2'])

    PlottingDoublePulses_Efield = libDatabase.FilterDatabaseGreaterThan(PlottingDoublePulses_Efield, Efield2_compare*(1.-FieldTolerance), libPlotting.IndexTable['Efield2'])
    #print Header+"Dataset for the field: "+str(Efield2_value*1e-9)+' V/nm'
    #print PlottingDoublePulses_Efield
    # TODO: as fields are all slighly different, we should project on the neareast neighboring value, up to 1E8 V/m. 

    #print Header+"Available data for the specific field E2="+str(Efield2_value*1E-9)+" V/nm."
    #tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9 = libPlotting.CountTheNumberOfLaserParameters(PlottingDoublePulses_Efield)

    ## prints the current state of selection and stores the values of sweeped parameters
    runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingDoublePulses_Efield)

    
    # Replaces wavelengths with photon energy
    if(EnergyAxis):  #very dirty
       wavelengths2_Show = np.divide(2*np.pi*c*hbar/e, wavelengths2) #eV #E=hbar omega = 2 pi c hbar / wavelength
    else:
       wavelengths2_Show = wavelengths2
    
    ## Sorting data before plotting. 
    # If using Energies, then we sort by energy. 
    # If using wavelengths, sorting is made by wavelengths. 
    lists = sorted(zip(*[wavelengths2_Show, Ne_final, pulses1, ElectronEnergy_max]))
    wavelengths2_sorted, Ne_final_sorted, pulses1_sorted, Energy_max_sorted = list(zip(*lists))
    
    if(not EnergyAxis):
        wavelengths2_sorted_show = np.multiply(1e9, wavelengths2_sorted)
    else: #if our wavelength is actually an energy
        wavelengths2_sorted_show = wavelengths2_sorted
    wpi_sorted = np.divide(Ne_final_sorted, np.multiply(1e0, pulses1_sorted)) #NOTE: this is limited to identical pulse durations!
    
    # This loop is run len(Efield2_uniq) times
    #for Efield2_compare in FieldOfInterest: #{{{
    logger.warning((Header+"Strange there are 1 FieldOfInterest, but we are getting once more in this loop. Efield2_compare: "+str(Efield2_compare*1E-9)))
    #LocalCondition = False
    #if(Efield2_compare==0):
        #LocalCondition = (Efield2_value == Efield2_compare)
    #else:
        #LocalCondition = (abs(Efield2_value - Efield2_compare)/Efield2_value <= FieldTolerance) #Many loops should not induce one new plot each time. 
    #if(LocalCondition): #{{{ #If Efield2_value is close enough to a field we would like to see, then we plot one new line. If we already had the value accepted once, we plot again?? No... So we should construct a set of data by merging all these values. 
    logger.info("Available data for the specific field E2="+str(Efield2_compare*1E-9)+" V/nm.")
    tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9 = libPlotting.CountTheNumberOfLaserParameters(PlottingDoublePulses_Efield)

    #### PLOTTING
    
    #if(not EnergyAxis):
    #plot_field, = ax1.plot(wavelengths2_sorted_show, wpi_sorted,'-^', antialiased=antialias, color=Colors[Efield2_compare], linewidth=LineWidthCoeff_abs*LineWidths[Efield2_compare], label=r'$\widetilde{w}_{\mathrm{PI}}$, I='+str(libKeldysh.FieldToIntensity(Efield2_value))+r' W/m$^2$') #BACKUP
    #plot_field, = ax1.plot(wavelengths2_sorted_show, wpi_sorted,'-^', antialiased=antialias, label=r'$\widetilde{w}_{\mathrm{PI}}$, I= ? W/m$^2$')
    plot_field, = ax1.plot(wavelengths2_sorted_show, wpi_sorted,'-^', label=r'$\widetilde{w}_{\mathrm{PI}}$, $\lambda_1$='+str(wavelength1*1E9)+' nm, $E_1$='+str(Field1*1E-9)+', $E_2$='+str(Efield2_compare*1E-9)+r' V/nm')
    #linewidth=LineWidthCoeff_abs*LineWidths[Efield2_compare] # BUG 
    plot1 += [plot_field]
    
    #plot_field, = ax1.plot(wavelengths1_Keldysh_Show, np.multiply(
            ##PeakWpi_Gruzdev, np.divide(np.add(n0, np.multiply(-1,PeakWpi_Gruzdev*pulses1_sorted[0])), n0)),
                #'--.'+Colors[Efield2_compare], label=r'$w_{\mathrm{PI}}\frac{n_0-n_{\mathrm{exc}}{n_0}$, I='+str(libKeldysh.FieldToIntensity(Efield2_value))+r' W/m$^2$')
    #plot1 += [plot_field]
    
    #}}} #end if
        
    del PlottingDoublePulses_Efield
  #}}} #end for
#}}} #end for
  
  # We separate TDDFT and Keldysh predictions now. 
  #ax2 = plt.subplot(plotB)
  
  #for Efield2_compare in FieldOfInterest: #{{{
    #print Header+"Building comparison with Keldysh-Gruzdev..."
    #order = 50 
    ##wavelengths1_min = np.min(wavelengths1); wavelengths1_max = np.max(wavelengths1)
    #wavelengths1_min = 200e-9; wavelengths1_max = 2100e-9
    #wavelengths1_Keldysh = np.linspace(wavelengths1_min, wavelengths1_max, 500)
    #PeakWpi_Gruzdev = libKeldysh.GenerateKeldyshDatabase(Egap_SI, EffectiveMass[Efield2_compare], wavelengths1_Keldysh, Efield2_compare, order)
    ##PeakWpi_Gulley = GenerateKeldyshGulleyDatabase(Egap_SI, EffectiveMass[Efield2_compare], wavelengths1_Keldysh, Efield2_value, order, 1.)
    ##PeakWpi_Gruzdev_Indirect = libKeldysh.GenerateKeldyshDatabase(Egap_SI_indirect, meff, wavelengths1_Keldysh, Efield2_value, order)
    
    #PeakWpi_Tunneling_Direct = np.clip(libKeldysh.KeldyshTunnelingLimit(Egap_SI, meff_Tunnel, wavelengths1_Keldysh, Efield2_compare, False), 1E36, 1E46)
    ##PeakWpi_Tunneling_Indirect = np.clip(libKeldysh.KeldyshTunnelingLimit(Egap_SI_indirect, meff_Tunnel, wavelengths1_Keldysh, Efield2_value), 1E36, 1E46)
    ##PeakWpi_Gruzdev=
    #PeakWpi_Effective = np.power(np.power(PeakWpi_Gruzdev, -1) + np.power(PeakWpi_Tunneling_Direct,-1),-1)
    ## Dirty: choosing energy or wavelengths. Sorting the results with this choice. 
    #if(EnergyAxis):
        #wavelengths1_Keldysh_Show = np.divide(2.*np.pi*c*hbar/e, wavelengths1_Keldysh) #eV
    #else:
        #wavelengths1_Keldysh_Show = np.multiply(1E9, wavelengths1_Keldysh)
    
    ## resorting the results to plot as function of energy if necessary
    #lists = sorted(itertools.izip(*[wavelengths1_Keldysh_Show, PeakWpi_Gruzdev]))
    #wavelengths1_Keldysh_Show, PeakWpi_Gruzdev_Show = list(itertools.izip(*lists)) 
    
    #plot_field, = ax2.plot(wavelengths1_Keldysh_Show, PeakWpi_Gruzdev_Show, '-', antialiased=antialias, color=Colors[Efield2_compare], linewidth=LineWidthCoeff_rel*LineWidthCoeff_abs*LineWidths[Efield2_compare], label=r'$\widetilde{w}_{\mathrm{PI}}$ [Gruzdev], I='+str(libKeldysh.FieldToIntensity(Efield2_value))+r' W/m$^2$')
    
    ##plot_field, = ax1.plot(wavelengths1_Keldysh_Show, PeakWpi_Gruzdev_Show, '-', linewidth=LineWidthCoeff_rel*LineWidthCoeff_abs*LineWidths[Efield2_compare], label=r'$\widetilde{w}_{\mathrm{PI}}$ [Gruzdev], I='+str(libKeldysh.FieldToIntensity(Efield2_value))+r' W/m$^2$')
    
    #plot1 += [plot_field]
  ##}}} #end for
    
  #if(EnergyAxis): 
    #ax1.set_xscale('log')
  ax1.set_yscale('log')
  ax1.set_ylabel(r'$\tilde{w}_{PI}$ (m$^{-3}$.s$^{-1}$)')
  #ax2.set_yscale('log')
  
  if(ShowDensities): #{{{
    ax12 = ax1.twinx()
    #plt.title('Final density as function of laser wavelength')
    logger.info("Total number of datasets: "+str(len(set(Efield2_uniq))))
    
    for Efield2_value in Efield2_uniq: #for each value of Efield1, we make a new search, and plot one new line in the chart
        #PlottingDoublePulses_Efield = libDatabase.FilterDatabase(PlottingDoublePulses, Efield2_value, IndexTable['Efield1']) #select the local Efield
        # TODO: how to merge data with values close to each others? Maybe introduce tolerance of 10% ?

        PlottingDoublePulses_Efield = libDatabase.FilterDatabaseLowerThan(PlottingDoublePulses, Efield2_value*(1.+FieldTolerance), libPlotting.IndexTable['Efield2'])

        PlottingDoublePulses_Efield = libDatabase.FilterDatabaseGreaterThan(PlottingDoublePulses_Efield, Efield2_value*(1.-FieldTolerance), libPlotting.IndexTable['Efield2'])
        
        #print PlottingDoublePulses_Efield
        ## prints the current state of selection and stores the values of sweeped parameters



        runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingDoublePulses_Efield)

        
        if(EnergyAxis): #very dirty!
          wavelengths2=np.divide(2.*np.pi*c*hbar/e, wavelengths2)
        
        lists = sorted(zip(*[wavelengths2, Ne_final, ElectronEnergy_max])) #NOTE: Ne_max makes more sense than Ne_final
        wavelengths2_sorted, Ne_final_sorted, Energy_max_sorted = list(zip(*lists))
        if(not EnergyAxis): 
            wavelengths2_sorted = np.multiply(1e9, wavelengths2_sorted)
        
        ax12.set_ylabel(r'Final density $n_{\mathrm{exc}}$ (m$^{-3}$)')
        #print Efield2_value
        # TODO: PUT here a select case which attributes a new color to a new field interval
        # OR USE A DICTIONNARY. 
        
        # SAME CODE, just more clever. 
        for Efield2_compare in FieldOfInterest:
            if(abs(Efield2_value - Efield2_compare)/Efield2_value < FieldTolerance):
                #if(not EnergyAxis):
                plot_field, = ax12.plot(wavelengths2_sorted, Ne_final_sorted, '^', color=Colors[Efield2_compare], label=r'$E_{peak}=$'+str(round(Efield2_compare*1E-9))+' V/nm')

                #plot_field, = ax12.plot(wavelengths2_sorted, Ne_final_sorted, '^', label=r'$E_{peak}=$'+str(round(Efield2_compare*1E-9))+' V/nm')
                plot1 += [plot_field]
                #, I='+str(libKeldysh.FieldToIntensity(Efield2_value))+r' W/m$^2$')
                break
        
        #if(1e8*(1.- FieldTolerance) < Efield2_value < 1e8*(1.+ FieldTolerance)):
            ##print Header+"Dataset for the field: "+str(Efield2_value*1e-9)+' V/nm'
            #plot_field, = ax12.plot(wavelengths2_sorted, Ne_final_sorted, '^r', label=r'$E_{peak}=$'+str(round(Efield2_value*1E-9))+' V/nm')
            #plot1 += [plot_field]
            ##, I='+str(libKeldysh.FieldToIntensity(Efield2_value))+r' W/m$^2$')
            #break
        #if(1e9*(1.- FieldTolerance) < Efield2_value < 1e9*(1.+ FieldTolerance)):
            ##print Header+"Dataset for the field: "+str(Efield2_value*1e-9)+' V/nm'
            #plot_field, = ax12.plot(wavelengths2_sorted, Ne_final_sorted, '^b', label=r'$E_{peak}=$'+str(round(Efield2_value*1E-9))+' V/nm') 
            #plot1 += [plot_field]
            ##, I='+str(libKeldysh.FieldToIntensity(Efield2_value))+r' W/m$^2$')
            #break
        #if(5e9*(1.- FieldTolerance) < Efield2_value < 5e9*(1.+ FieldTolerance)):
            ##print Header+"Dataset for the field: "+str(Efield2_value*1e-9)+' V/nm'
            #plot_field, = ax12.plot(wavelengths2_sorted, Ne_final_sorted, '^g', label=r'$E_{peak}=$'+str(round(Efield2_value*1E-9))+' V/nm')
            #plot1 += [plot_field]
            ##, I='+str(libKeldysh.FieldToIntensity(Efield2_value))+r' W/m$^2$')
            #break
        #if(1e10*(1.- FieldTolerance) < Efield2_value < 1e10*(1.+ FieldTolerance)):
            ##print Header+"Dataset for the field: "+str(Efield2_value*1e-9)+' V/nm'
            #plot_field, = ax12.plot(wavelengths2_sorted, Ne_final_sorted, '^k', label=r'$E_{peak}=$'+str(round(Efield2_value*1E-9))+' V/nm')
            #plot1 += [plot_field]
            ##, I='+str(libKeldysh.FieldToIntensity(Efield2_value))+r' W/m$^2$')
            #break
        del PlottingDoublePulses_Efield
  #}}}  
  labels1 = [l.get_label() for l in plot1]
  ax1.legend(plot1, labels1, loc='best') #outside right, middle
  #ax1.legend(loc=3, borderaxespad=0., bbox_to_anchor=(1.05, 0.05))
  if(not EnergyAxis): 
    ax1.set_xlabel('Laser wavelength (nm)')
    ax1.set_xlim([200.,2100.])
    #ax2.set_xlabel('Laser wavelength (nm)')
    #ax2.set_xlim([200.,2100.])
  else: 
    ax1.set_xlabel('Photon energy (eV)')
    ax1.set_xlim([0.5,4.0])
    #ax2.set_xlabel('Photon energy (eV)')
    #ax2.set_xlim([0.5,4.0])
  ax1.set_ylim([1E35,1E43])
  #ax2.set_ylim([1E35,1E43])
  #ax1.set_xscale('log')
  #ax1.set_yscale('log')
  #ax1.grid()
  
  ResonantWavelength = h*c/Egap_SI
  logger.info("Adding vertical bar for the resonant transitions at: "+str(round(1e9*ResonantWavelength))+"nm.")
  
  if(not EnergyAxis):
    ax1.axvline(x=1e9*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax1.axvline(x=1e9*2*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax1.axvline(x=1e9*3*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax1.axvline(x=1e9*4*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax1.axvline(x=1e9*wavelength1, color='k', linestyle='--', linewidth=1)
    
    #ax2.axvline(x=1e9*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    #ax2.axvline(x=1e9*2*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    #ax2.axvline(x=1e9*3*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    #ax2.axvline(x=1e9*4*ResonantWavelength, color='k', linestyle='--', linewidth=1)
  else:
    ax1.axvline(x=Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax1.axvline(x=1./2.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax1.axvline(x=1./3.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax1.axvline(x=1./4.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax1.axvline(x=2.*np.pi*c/wavelength1/e, color='r', linestyle='-', linewidth=2)
    
    #ax2.axvline(x=Egap_SI/e, color='k', linestyle='--', linewidth=1)
    #ax2.axvline(x=1./2.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    #ax2.axvline(x=1./3.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    #ax2.axvline(x=1./4.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
  
  ## Saving this first double figure
  #plt.tight_layout()
  filename = 'Metaplot-DoublePulse-BiColor-EffectOfWavelength'+str(round(wavelength1*1E9,0))+'nm-Efield1_Tol-'+str(round(Field1*1E-9,1))+'Vnm-Keldysh-'+str(int(pulseDuration*1E15))+"-CEP"+str(round(CEP2_request*180./np.pi,0))
  plt.savefig(filename+'.eps')
  plt.savefig(filename+'.png')
  
  ## Computes the free carrier energy contained in the volume <SampleVolume>.
  def EnergyFromFreeCarrierAbsorption(wavelengths1_sorted, pulseDuration, Efield, Ne_final_sorted, w_PI_sorted, Egap_SI = 2.56*e, meff = 0.2226, SampleVolume=2 * 5E28**-1): #{{{
    # dEnergy/dt = alpha * intensity + ( n har omega - Egap ) x sigma_n I^n
    # where we approach sigma_n I^n by pulse-averaged w_PI 
    # Energy ~ alpha * intensity + ( n har omega - Egap ) x w_PI * tau
    omega   = np.divide(2. * np.pi * c, wavelengths1_sorted)
    #print Header+"omega: "+str(omega) #Valid
    nu = 6e-15**-1 #1.1 fs: Solowski-Tinten adjustment. Can also use the one from Shunsuke Sato, PRB 2015 (Numerical pump probe)
    logger.info("nu: "+str(nu)) #Valid
    #epsilon_Si0 = 13.64 + 0.048j #TODO: get the optical data in extended spp theory for Si
    epsilon_Si0 = libOpticalData.ImportPalikDatabase_epsilon_fromNK(wavelengths1_sorted) ##TODO: replace by LDA !
    logger.debug("Dielectric permittivity of Si at 800 nm : "+str(libOpticalData.ImportPalikDatabase_epsilon_fromNK(800e-9)))
    omega_p2    = np.multiply(Ne_final_sorted, e**2 / m_e / meff / epsilon_0)
    #print "omega_p: "+str(np.sqrt(omega_p2)) # Valid
    # NOTE: THIS IS A CHOICE. WE JUST HAVE TO BE CONSISTENT ON WHICH PROCESS WE WISH FOR. 
    Epsilon = 1E0 - np.divide(omega_p2, np.power(omega, 2)) * np.divide(1.,1. + np.divide(1.j * nu, omega))
    #Epsilon = epsilon_Si0 - np.divide(omega_p2, np.power(omega, 2)) * np.divide(1.,1. + np.divide(1.j * nu, omega))
    #NumberOfPhotons = Egap_SI / (hbar * omega)
    NumberOfPhotons = np.divide(Egap_SI/hbar, omega)
    logger.info(NumberOfPhotons)
    #NumberOfPhotons_round = np.round(NumberOfPhotons,0)
    NumberOfPhotons_round = np.ceil(NumberOfPhotons)
    logger.debug(Header+"sqrt(Epsilon): "+str(np.min(np.sqrt(Epsilon)))+", "+str(np.max(np.sqrt(Epsilon)))) #Valid
    alpha                = np.multiply(2.* omega/c, np.imag(np.sqrt(Epsilon))) #valid
    EnergyExcess         = np.multiply(NumberOfPhotons_round * hbar * omega, Ne_final_sorted) * SampleVolume / e  #result in eV
    EnergyBremsstrahlung = np.multiply( omega * epsilon_0 * (Epsilon.imag), np.multiply(Efield, Efield)) * SampleVolume / e  * pulseDuration
    EnergyBremsstrahlung_approx =  EnergyBremsstrahlung #NOTE: not used for now. Another collision freq? 
    #EnergyBremsstrahlung_approx = np.multiply( 4.*np.pi / c * np.divide(Epsilon.imag, (np.sqrt(Epsilon)).real), )) ) * SampleVolume / e
    logger.info(Header+"FcrEnergy (eV): "+str(np.min(EnergyBremsstrahlung))+", "+str(np.max(EnergyBremsstrahlung)))
    logger.info(Header+"Excess energy (n hbar omega) (eV): "+str(np.min(EnergyExcess))+", "+str(np.max(EnergyExcess)))
    EnergyFreeCarrier    = EnergyExcess + EnergyBremsstrahlung
    EnergyFreeCarrier_approx = EnergyExcess + EnergyBremsstrahlung_approx
    return EnergyFreeCarrier, EnergyFreeCarrier_approx
  #}}}
  
  ## TODO: TDDFT vs Drude should become another figure. We don't need to have them side by side. 
  del labels1; del plot1
  plt.figure(figsize=(12,6))
  
  plot1=[]
  
  ax3 = plt.subplot(111)
  #plt.title('Irreversibility factor')
  #plt.plot(libKeldysh.FieldToIntensity(Efield1), irreversibility, 'o', label='Irreversibility')
  for Efield2_compare in FieldOfInterest: #for each value of Efield1, we make a new search, and plot one new line in the chart
    #PlottingDoublePulses_Efield = libDatabase.FilterDatabase(PlottingDoublePulses, Efield2_value, IndexTable['Efield1']) #select the local Efield




    PlottingDoublePulses_Efield = libDatabase.FilterDatabaseLowerThan(PlottingDoublePulses, Efield2_compare*(1.+FieldTolerance), libPlotting.IndexTable['Efield2'])

    PlottingDoublePulses_Efield = libDatabase.FilterDatabaseGreaterThan(PlottingDoublePulses_Efield, Efield2_compare*(1.-FieldTolerance), libPlotting.IndexTable['Efield2'])
    #print Header+"Dataset for the field: "+str(Efield2_value*1e-9)+' V/nm'
    #print PlottingDoublePulses_Efield
    # TODO: as fields are all slighly different, we should project on the neareast neighboring value, up to 1E8 V/m. 
    
    ## prints the current state of selection and stores the values of sweeped parameters



    runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingDoublePulses_Efield)

    
    # 
    if(EnergyAxis):
        wavelengths2_Show = np.divide(2.*np.pi*c*hbar/e, wavelengths2) #eV
    else: 
        wavelengths2_Show = np.multiply(1e9, wavelengths2)
    
    lists = sorted(zip(*[wavelengths2_Show, wavelengths2, Ne_final, W_pi, ElectronEnergy_max])) #NOTE: Ne_max makes more sense than Ne_final
    wavelengths2_sorted_show, wavelengths2_sorted, Ne_final_sorted, W_pi_sorted, Energy_max_sorted = list(zip(*lists))
    #print Efield2_value
    
    #for Efield2_compare in FieldOfInterest: 
    #LocalCondition = False
    #if(Efield2_compare==0):
        #LocalCondition = (Efield2_value == Efield2_compare)
    #else:
        #LocalCondition = (abs(Efield2_value - Efield2_compare)/Efield2_compare <= FieldTolerance)
    #if(LocalCondition):  #if we have captured the fields we are interested in, we plot  
    #ax3.plot(wavelengths2_sorted_show, Energy_max_sorted, '-', linewidth=LineWidthCoeff_abs*LineWidths[Efield2_compare], color=Colors[Efield2_compare])
    plot_field, = ax3.plot(wavelengths2_sorted_show, Energy_max_sorted, '-^', label=r'$\lambda_1=$'+str(wavelength1*1E9)+' nm, $E_1=$'+str(Field1*1E-9)+r' V/nm, $E_2=$'+str(Efield2_compare*1E-9)+" V/nm.")
    logger.info((Header+"Plotting energy: Adding one plot..."))
    logger.info(wavelengths2_sorted_show)
    logger.debug(Energy_max_sorted)
    plot1 += [plot_field]
    #ax3.plot(wavelengths2_sorted_show, Energy_max_sorted, '-')
            
    del PlottingDoublePulses_Efield
    
  #plt.xlabel(r'Intensity (W/m$^2$)')
  if(not EnergyAxis):
    ax3.set_xlabel('Laser wavelength (nm)')
  else: 
    ax3.set_xlabel('Photon energy (eV)')
  #plt.ylabel('Irreversibility factor')
  ax3.set_ylabel('Energy absorbed by electrons (eV)')
  #plt.ylim((0,1))
  #ax2.set_xscale('log')
  
  ## TODO: this can be removed, since irreversibility is gauge-dependent. 
  #if(PlotIrreversibility): #{{{
    #ax22 = ax2.twinx()
    #for Efield2_value in Efield2_uniq: #for each value of Efield1, we make a new search, and plot one new line in the chart
        ##PlottingDoublePulses_Efield = libDatabase.FilterDatabase(PlottingDoublePulses, Efield2_value, IndexTable['Efield1']) #select the local Efield
        #PlottingDoublePulses_Efield = libDatabase.FilterDatabaseGreaterThan(PlottingDoublePulses_Efield, Efield2_value*(1.-FieldTolerance), IndexTable['Efield1'])
        ##print Header+"Dataset for the field: "+str(Efield2_value*1e-9)+' V/nm'
        ##print PlottingDoublePulses_Efield
        ## TODO: as fields are all slighly different, we should project on the neareast neighboring value, up to 1E8 V/m. 
        
        ### prints the current state of selection and stores the values of sweeped parameters
        #runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime, FinalTime = libPlotting.ExtractData(PlottingDoublePulses_Efield)
        
        #lists = sorted(itertools.izip(*[wavelengths1, Ne_final, Energy_max, irreversibility])) #NOTE: Ne_max makes more sense than Ne_final
        #wavelengths1_sorted, Ne_final_sorted, Energy_max_sorted, irreversibility_sorted = list(itertools.izip(*lists))
        #wavelengths1_sorted_show = np.multiply(1e9, wavelengths1_sorted)
        
        #for Efield_compare in FieldOfInterest:
            #if(abs(Efield2_value - Efield2_compare)/Efield2_compare < FieldTolerance):
                #ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+'+Colors[Efield2_compare], label='Irreversibility')
        
        ###print Efield2_value
        ### TODO: PUT here a select case which attributes a new color to a new field interval
        ###if(1e8*(1.-FieldTolerance) < Efield2_value < 1e8*(1.+FieldTolerance)):
            ###ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+r', label='Irreversibility')
            ###ax22.plot(wavelengths1_sorted_show, libKeldysh.gammaKeldysh(Egap_SI, 0.228, Efield2_value, wavelengths1_sorted), '--r', label='Adiabaticity')
        ##if(1e9*(1.-FieldTolerance) < Efield2_value < 1e9*(1.+FieldTolerance)):
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+b', label='Irreversibility')
            ###ax22.plot(wavelengths1_sorted_show, libKeldysh.gammaKeldysh(Egap_SI, 0.228, Efield2_value, wavelengths1_sorted), '--b', label='Adiabaticity')
        ##if(5e9*(1.-FieldTolerance) < Efield2_value < 5e9*(1.+FieldTolerance)):
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+g', label='Irreversibility')
            ###ax22.plot(wavelengths1_sorted_show, libKeldysh.gammaKeldysh(Egap_SI, 0.228, Efield2_value, wavelengths1_sorted), '--g', label='Adiabaticity')
        ##if(1e10*(1.-FieldTolerance) < Efield2_value < 1e10*(1.+FieldTolerance)):
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+k', label='Irreversibility')
            ###ax22.plot(wavelengths1_sorted_show, libKeldysh.gammaKeldysh(Egap_SI, 0.228, Efield2_value, wavelengths1_sorted), '--k', label='Adiabaticity')
        #del PlottingDoublePulses_Efield
    
    
    #plt.ylabel(r'Irreversibility $\zeta$')
    #ax22.set_yscale('log')
  ##}}}
  
  # This part adds an analytical estimation based on Drude(n_exc_TDDFT)
  if(PlotFreeCarrierAbsorption):
      logger.info(Header+"Now plotting the free carrier energy absorption...")
      for Efield2_value in Efield2_uniq: #for each value of Efield1, we make a new search, and plot one new line in the chart
        #PlottingDoublePulses_Efield = libDatabase.FilterDatabase(PlottingDoublePulses, Efield2_value, libPlotting.IndexTable['Efield1']) #select the local Efield
        PlottingDoublePulses_Efield = libDatabase.FilterDatabaseLowerThan(PlottingDoublePulses, Efield2_value*(1.+FieldTolerance), libPlotting.IndexTable['Efield2'])

        PlottingDoublePulses_Efield = libDatabase.FilterDatabaseGreaterThan(PlottingDoublePulses_Efield, Efield2_value*(1.-FieldTolerance), libPlotting.IndexTable['Efield2'])
        #print Header+"Dataset for the field: "+str(Efield2_value*1e-9)+' V/nm'
        #print PlottingDoublePulses_Efield
        # TODO: as fields are all slighly different, we should project on the neareast neighboring value, up to 1E8 V/m. 
        
        ## prints the current state of selection and stores the values of sweeped parameters
        runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingDoublePulses_Efield)

        
        if(EnergyAxis):
            wavelengths2_Show = np.divide(2.*np.pi*c*hbar/e, wavelengths2) #eV
        else: 
            wavelengths2_Show = np.multiply(1e9, wavelengths2)
    
        lists = sorted(zip(*[wavelengths2_Show, wavelengths2, Ne_final, W_pi, ElectronEnergy_max])) #NOTE: Ne_max makes less sense than Ne_final
        wavelengths2_sorted_show, wavelengths2_sorted, Ne_final_sorted, W_pi_sorted, Energy_max_sorted = list(zip(*lists))
    #print Efield2_value
    
        
        #lists = sorted(itertools.izip(*[wavelengths2, Ne_final, W_pi, Energy_max, irreversibility])) #NOTE: Ne_max makes less sense than Ne_final
        #wavelengths2_sorted, Ne_final_sorted, W_pi_sorted, Energy_max_sorted, irreversibility_sorted = list(itertools.izip(*lists))
        #wavelengths2_sorted_show = np.multiply(1e9, wavelengths2_sorted)
        
        for Efield2_compare in FieldOfInterest:
            LocalCondition = False
            if(Efield2_compare==0):
                LocalCondition = (Efield2_value == Efield2_compare)
            else:
                LocalCondition = (abs(Efield2_value - Efield2_compare)/Efield2_compare <= FieldTolerance)
            if(LocalCondition):
                EnergyAnalytical, EnergyAnalytical_approx = EnergyFromFreeCarrierAbsorption(pulses1, wavelengths2_sorted, Efield2_value, Ne_final_sorted, W_pi_sorted, Egap_SI, EffectiveMass[Efield2_compare])
                                
                plot_field, = ax3.plot(wavelengths2_sorted_show, EnergyAnalytical, '--', color=Colors[Efield2_compare], label=r'$E_1=$'+str(Field1*1E-9)+r' V/nm, $E_2=$'+str(Efield2_compare*1E-9)+" V/nm.")
                #ax2.plot(wavelengths2_sorted_show, EnergyAnalytical, '--', label='Free carrier absorption')
                plot1 += [plot_field]
            
        ##if(1e8*(1.-FieldTolerance) < Efield2_value < 1e8*(1.+FieldTolerance)):
            ##print Header+"Found field 0.1 V/nm"
            ##EnergyAnalytical, EnergyAnalytical_approx = EnergyFromFreeCarrierAbsorption(pulses1, wavelengths1_sorted, Efield2_value, Ne_final_sorted, W_pi_sorted, Egap_SI)
            ##ax2.plot(wavelengths1_sorted_show, EnergyAnalytical, '--+r', label='Free carrier absorption')
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+r', label='Irreversibility')
            ##ax22.plot(wavelengths1_sorted_show, libKeldysh.gammaKeldysh(Egap_SI, 0.228, Efield2_value, wavelengths1_sorted), '--r', label='Adiabaticity')
        #if(1e9*(1.-FieldTolerance) < Efield2_value < 1e9*(1.+FieldTolerance)):
            #print Header+"Found field 1.0 V/nm"
            #EnergyAnalytical, EnergyAnalytical_approx = EnergyFromFreeCarrierAbsorption(pulses1, wavelengths1_sorted, Efield2_value, Ne_final_sorted, W_pi_sorted, Egap_SI)
            #ax2.plot(wavelengths1_sorted_show, EnergyAnalytical, '--+b', label='Free carrier absorption')
            ##ax2.plot(wavelengths1_sorted_show, EnergyAnalytical_approx, '--b', label='Free carrier absorption')
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+b', label='Irreversibility')
            ##ax22.plot(wavelengths1_sorted_show, libKeldysh.gammaKeldysh(Egap_SI, 0.228, Efield2_value, wavelengths1_sorted), '--b', label='Adiabaticity')
        #if(5e9*(1.-FieldTolerance) < Efield2_value < 5e9*(1.+FieldTolerance)):
            #print Header+"Found field 5 V/nm"
            #EnergyAnalytical, EnergyAnalytical_approx = EnergyFromFreeCarrierAbsorption(pulses1, wavelengths1_sorted, Efield2_value, Ne_final_sorted, W_pi_sorted, Egap_SI)
            #ax2.plot(wavelengths1_sorted_show, EnergyAnalytical, '--+g', label='Free carrier absorption')
            ##ax2.plot(wavelengths1_sorted_show, EnergyAnalytical_approx, '--g', label='Free carrier absorption')
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+g', label='Irreversibility')
            ##ax22.plot(wavelengths1_sorted_show, libKeldysh.gammaKeldysh(Egap_SI, 0.228, Efield2_value, wavelengths1_sorted), '--g', label='Adiabaticity')
        #if(1e10*(1.-FieldTolerance) < Efield2_value < 1e10*(1.+FieldTolerance)):
            #print Header+"Found field 10 V/nm"
            #EnergyAnalytical, EnergyAnalytical_approx = EnergyFromFreeCarrierAbsorption(pulses1, wavelengths1_sorted, Efield2_value, Ne_final_sorted, W_pi_sorted, Egap_SI)
            #ax2.plot(wavelengths1_sorted_show, EnergyAnalytical, '--+k', label='Free carrier absorption')
            ##ax2.plot(wavelengths1_sorted_show, EnergyAnalytical_approx, '--k', label='Free carrier absorption')
            ##ax22.plot(wavelengths1_sorted_show, irreversibility_sorted, '--+k', label='Irreversibility')
            ##ax22.plot(wavelengths1_sorted_show, libKeldysh.gammaKeldysh(Egap_SI, 0.228, Efield2_value, wavelengths1_sorted), '--k', label='Adiabaticity')
        del PlottingDoublePulses_Efield
    
    #plt.ylabel(r'Irreversibility $\zeta$')
    #ax22.set_yscale('log')
      #ax22 = ax2.twinx()
  #plt.ylabel(r'Adiabaticity $\gamma$, Irreversibility $\zeta$')
  #ax22.set_ylim((0,1))
  labels1 = [l.get_label() for l in plot1]
  ax3.legend(plot1, labels1, loc='best') #outside right, middle
  #ax3.set_xscale("log")
  ax3.set_yscale("symlog", linthreshy=1e-9)
  #ax2.set_ylim((1E-9,5E2))
  #ax22.grid()
  if(not EnergyAxis): 
    plt.xlim((200.,2100.))
  else: 
    plt.xlim((0.5,7.))
  
  
  #ResonantWavelength = h*c/Egap_SI
  logger.info("** Adding vertical bar for the resonant transitions at: "+str(round(1e9*ResonantWavelength))+"nm.")
  
  if(not EnergyAxis):
    ax3.axvline(x=1e9*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax3.axvline(x=1e9*2*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax3.axvline(x=1e9*3*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax3.axvline(x=1e9*4*ResonantWavelength, color='k', linestyle='--', linewidth=1)
    ax3.axvline(x=1E9*wavelength1, color='r', linestyle='-', linewidth=2)
  else:
    ax3.axvline(x=Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax3.axvline(x=1./2.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax3.axvline(x=1./3.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax3.axvline(x=1./4.*Egap_SI/e, color='k', linestyle='--', linewidth=1)
    ax3.axvline(x=2.*np.pi*c/wavelength1/e, color='r', linestyle='-', linewidth=2)
  #ax2.legend(loc=2, borderaxespad=0.)
  #ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
  #ax22.legend(loc=1, borderaxespad=0.)

  # Saving the plot
  #plt.tight_layout()
  #filename = 'Metaplot-DoublePulse-BiColor-EffectOfWavelength'+str(round(wavelength1*1E9,0))+'nm-Efield1_Tol-'+str(round(Efield1*1E-9,0))+'Vnm-Keldysh-'+str(int(pulseDuration*1E15))
  filename='Metaplot-DoublePulse-BiColor-EffectOfWavelength'+str(round(wavelength1*1E9,0))+'nm-Efield1_Tol-'+str(round(Field1*1E-9,1))+'Vnm-Drude-'+str(int(pulseDuration*1E15))
  plt.savefig(filename+'.eps')
  plt.savefig(filename+'.png')
  logger.info((Header+"Bicolor_plotEffectOfWavelength_Tolerance: successfully wrote the bicolor Stark analysis."))
  if(ShowFinalPlot==1):
    plt.show()

  del ax1, PlottingDoublePulses
  #}}}

## Plots Nexc | Energy (CEP) for bicolor simulations. Selection is made with function parameters. A 1 percent tolerance is used on the <Efield1> and <Efield2> amplitudes. The rest is exactly tested. 
def Bicolor_plotEffectOfCEP(PlottingTable, Color1=800e-9, duration1=20e-15, FieldAmplitude1=1E9, Color2=1600e-9, duration2=20e-15, FieldAmplitude2=1E9, PulseShape=3, pulseRiseTime=0, ShowFinalPlot=0, Egap_SI=2.56*e):
  print("#====== Plot: DOUBLE COLOR PULSES, EFFECT OF THE CEP, SAME POLARIZATIONS, SAME AMPLITUDE: final density as function of field ======")
  #Filter: double pulses, same wavelengths, same polarization
  
  # Simple pulse, single color
  PlottingSinglePulses = PlottingTable #get the full database
  PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, libPlotting.IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
  PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, libPlotting.IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select double pulses.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0., libPlotting.IndexTable['CrossPolarization']) #removes fields where CrossPolarization is greater than 0.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, Color1, libPlotting.IndexTable['wavelengths1']) #keeps only 800 nm wavelenghts. #TODO: this field should be making a loop in available wavelengths given in wavelengths1_uniq.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, Color2, libPlotting.IndexTable['wavelengths2']) #keeps only 800 nm wavelenghts. #TODO: this field should be making a loop in available wavelengths given in wavelengths1_uniq.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0, libPlotting.IndexTable['FreezeHXC']) #keeps only non-frozen XC.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, PulseShape, libPlotting.IndexTable['PulseType']) #keeps only non-frozen XC.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, duration1, libPlotting.IndexTable['pulses1']) #removes fields where tau != pulseDuration
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, duration2, libPlotting.IndexTable['pulses2']) #removes fields where tau != pulseDuration
  PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, FieldAmplitude1*1.01, libPlotting.IndexTable['Efield1'])
  PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, FieldAmplitude1*0.99, libPlotting.IndexTable['Efield1'])
  PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, FieldAmplitude2*1.01, libPlotting.IndexTable['Efield2'])

  PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, FieldAmplitude2*0.99, libPlotting.IndexTable['Efield2'])
  
  if(pulseRiseTime != 0): 
    PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, pulseRiseTime, libPlotting.IndexTable['RiseTime']) #removes fields where tau != pulseDuration

  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = libPlotting.CountTheNumberOfLaserParameters(PlottingSinglePulses)
  
  plt.figure(figsize=(8,8./GoldNumber))
  ax1 = plt.subplot(111) #111 or plotA if using plotB
  plt.title(r'$n_{\mathrm{exc}}(\phi_{\mathrm{CEP}})$')
  # Number of sweeped parameters: 3 (Efield, Delay, CrossPolarization)
  #plt.scatter(libKeldysh.FieldToIntensity(Efield1), Ne_final, c=Delay, s=(30.+CrossPolarization*50.), label=r'$N_{exc}$[Field]')
  plot1 = [] #np.empty((0,0))
  #for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
  #PlottingSinglePulses_CEP = libDatabase.FilterDatabase(PlottingSinglePulses, CEP1_value, IndexTable['CEP1']) #select the local CEP
  
  ## prints the current state of selection and stores the values of sweeped parameters
  runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max,Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingSinglePulses)

  
  #plot_field, = ax1.plot(libKeldysh.FieldToIntensity(Efield1), Ne_final, '-', label=r'$N_{exc}^{final}$, $\phi_{CE}$='+str(round(CEP1_value*180./pi)))
  #lists = sorted(itertools.izip(*[libKeldysh.FieldToIntensity(Efield1), Ne_final]))
  lists = sorted(zip(*[CEP2, Ne_final, W_pi, ElectronEnergy_max, Efield2]))
  CEP2_sorted, Ne_final_sorted, W_pi_sorted, Energy_max_sorted, Efield2_sorted = list(zip(*lists))
  plot_field, = ax1.plot(np.multiply(CEP2_sorted,180./pi), Ne_final_sorted, 'r-+', label=r'$n_{\mathrm{exc}}^{\mathrm{final}}$, $E^{\mathrm{peak}}_{1,2}=$'+str(round(FieldAmplitude1*1E-9))+' V/nm.')
  logger.debug(plot_field)
  logger.debug(Header+"** Debug info:"+str(np.shape(plot_field)))
  plot1 += [plot_field] # does cumulate each CEP now. Type: Line2D, not np.array. 
    
  #ax1.legend(loc=4, borderaxespad=0.) #bbox_to_anchor=(1.05, 1), 
  #plt.xlabel(r'Peak intensity (W/m$^2$)')
  plt.xlabel(r'Carrier envelope phase $\phi$ (deg)')
  plt.ylabel(r'Final density (m$^{-3}$)')
  #  plt.ylim((1E16,4.*5e28))
  #ax1.set_xscale('log')
  #ax1.set_yscale('log')
  #plot1_list = plot1.tolist()
    #print Efield1, CEP1_value
    
  #Adding axis of gammaKeldysh, but for all CEPs !!
  #runid, wavelengths1, pulses1, Efield1, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, Functional, FreezeHXC, Ne_max, Ne_final, irreversibility, TotalCurrent_final, Energy_max, W_pi, PulseType, RiseTime = libPlotting.ExtractData(PlottingSinglePulses)
  
  ax12 = ax1.twinx()
  #for CEP1_value in CEP1_uniq: #for each value of CEP1, we make a new search, and plot one new line in the chart
  #PlottingSinglePulses_CEP = libDatabase.FilterDatabase(PlottingSinglePulses, CEP1_value, IndexTable['CEP1']) #select the local CEP
  
  ## prints the current state of selection and stores the values of sweeped parameters



  runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max,  ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max,Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingSinglePulses)

  
  #plot_field, = ax1.plot(libKeldysh.FieldToIntensity(Efield1), Ne_final, '-', label=r'$N_{exc}^{final}$, $\phi_{CE}$='+str(round(CEP1_value*180./pi)))
  #lists = sorted(itertools.izip(*[libKeldysh.FieldToIntensity(Efield1), Ne_final]))
  lists = sorted(zip(*[CEP2, Ne_final, W_pi, ElectronEnergy_max, pulses2, Efield2]))
  CEP2_sorted, Ne_final_sorted, W_pi_sorted, Energy_max_sorted, pulses1_sorted, Efield2_sorted  = list(zip(*lists))
  #plt.plot(np.multiply(CEP2_sorted,180./pi), np.divide(Ne_final_sorted, pulses1_sorted), '-') #, label=r'$N_{exc}/\tau$') #, I='+str(libKeldysh.FieldToIntensity(Efield1_value))+r' W/m$^2$')
  plot_field, = ax12.plot(np.multiply(CEP2_sorted,180./pi), Energy_max_sorted, 'b-', label=r'$\xi_{\mathrm{el}}^{\mathrm{TDDFT}}$ , $E^{\mathrm{peak}}_{1,2}=$'+str(round(FieldAmplitude2*1E-9,0))+' V/nm')
  plot1 += [plot_field]
  
  #ax12.set_yscale('log')
  ax12.set_ylabel(r'$\xi_{\mathrm{el}}^{\mathrm{TDDFT}}$ (eV)')
  #ax.ylabel('Virtual exc. factor')
  
  labels1 = [l.get_label() for l in plot1]
  ax1.legend(plot1, labels1, loc='best', fancybox=True, framealpha=0)
  
  
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
  ExportName = "EffectOfCEP-E1E2-"+str(round(FieldAmplitude1*1E-9,0))+'-tau1tau2-'+str(round(duration1*1E15,0))
  plt.savefig('Metaplot-Bicolor-'+ExportName+'.eps')
  plt.savefig('Metaplot-Bicolor-'+ExportName+'.png')
  if(ShowFinalPlot==1):
    plt.show()

  del ax1, PlottingSinglePulses, wavelength_uniq
#}}}

## Plots the N_exc(wavelength, Efield) (2D plot)
# Not usable, unless you have computed a continuum of wavelengths AND intensities (crazy amount of computer time - 1M core-hour per figure)
# VALID, but of course interference of 2 waves is oscillating very much. 
def plotEffectOfWavelength_2dMap(PlottingTable, pulseDuration=20e-15, PulseType=2, RiseCycles=0, ShowFinalPlot=0, Egap_SI=2.56*e):
  print("#====== Plot: SINGLE PULSE, EFFECT OF THE WAVELENGTH, SAME CEP, SAME POLARIZATIONS: final density as function of field ======")
  #Filter: double pulses, same wavelengths, same polarization

  # Simple pulse, single color
  PlottingSinglePulses = PlottingTable #get the full database
  PlottingSinglePulses = libDatabase.FilterDatabaseGreaterThan(PlottingSinglePulses, 1E0, libPlotting.IndexTable['Efield1']) #removes fields where Epeak1 < 1 V/m. #works
  PlottingSinglePulses = libDatabase.FilterDatabaseLowerThan(PlottingSinglePulses, 1E0, libPlotting.IndexTable['Efield2']) #removes fields where Epeak2 > 1 V/m to select single pulses.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0., libPlotting.IndexTable['CrossPolarization']) #removes fields where CrossPolarization is different from 0.
  PlottingSinglePulses = libDatabase.FilterDatabase(PlottingSinglePulses, 0., libPlotting.IndexTable['CEP1']) #keeps only CEP = 0.


  ## TODO: check if CEP could count at smaller wavelength. 

  logger.info(Header+"Available data: ")
  wavelength_uniq, CEP1_uniq, TheoryLevel_uniq, Efield1_uniq, pulses_uniq,wavelength2_uniq,CEP2_uniq,Efield2_uniq,pulses2_uniq = libPlotting.CountTheNumberOfLaserParameters(PlottingSinglePulses)

  
  #plt.figure(figsize=(SizeX,SizeY)) #new figure
  DensityMax2D = np.empty((0,3)) #preparing the 2D plot
  #print DensityMax2D
  #ax1 = plt.subplot(111)
  plt.title('Final density as function of laser wavelength')
  logger.info(Header+"Total number of datasets: "+str(len(set(Efield1_uniq))))
  
  ##Building up the stack to plot here, therefore we could plot anytype of plot. 
  for Efield1_value in Efield1_uniq: #for each value of Efield1, we make a new search, and plot one new line in the chart
    PlottingSinglePulses_Efield = libDatabase.FilterDatabase(PlottingSinglePulses, Efield1_value, libPlotting.IndexTable['Efield1']) #select the local Efield
    for wavelength_value in wavelength_uniq:
      PlottingSinglePulses_Efield_wavelength = libDatabase.FilterDatabase(PlottingSinglePulses_Efield, wavelength_value, libPlotting.IndexTable['wavelengths1']) #select the local Efield
      runid, Machine, RunDate, kpointgrid, realmeshgrid, TimeStep, wavelengths1, pulses1, Efield1, polarization_angle, CEP1, wavelengths2, pulses2, Efield2, CEP2, Delay, CrossPolarization, TheoryLevel, FieldInduced, GaugeFieldAmplitude, GaugeFieldDelay, Functional, FreezeHXC, MoveIon, Ne_max, Ne_final, irreversibility, TotalCurrent_final, TotalCurrent_max, TotalCurrent_AU_Xmax, TotalCurrent_AU_Ymax, TotalCurrent_AU_Zmax, ElectronEnergy_max, EigenvaluesEnergy_max, HartreeEnergy_max, ExchangeEnergy_max, ElectronEnergy_afterpulse, EigenvaluesEnergy_afterpulse, HartreeEnergy_afterpulse, ExchangeEnergy_afterpulse, EnergyPoynting_SI_sum_final, Eext_SI_max, Ne_afterpulse, W_pi, Wpi_afterpulse, PulseType, RiseTime, FinalTime, timeToFilterFrom , LastRunTimeStepDuration , NumberOfProcesses = libPlotting.ExtractData(PlottingSinglePulses_Efield_wavelength)

      
      if(len(Ne_final) > 0e0): #Filtering system is weird
        #DensityMax2D=np.vstack((DensityMax2D, [wavelengths1[0], Efield1[0], Ne_final[0]]))
        ToBeAdded=np.array([wavelengths1[0], Efield1[0], Ne_final[0]])
        #print "Stacking..."
        #print ToBeAdded
        DensityMax2D=np.vstack((DensityMax2D, ToBeAdded))
  
  # Stack is built. Export, show and plot on a heat map. 
  logger.info("** Exporting the sub-plot2D: Ne_final(wavelength, Field)... Useful to visualize available data with GNUplot. ")
  libDatabase.ExportToTxt(DensityMax2D, "Wavelength-Field-NeFinal.csv")

  plt.figure()
  #CS = plt.imshow(DensityMax2D, cmap=plt.cm.Blues) #WORKS, but data are too small (no interpolation)
  CS = plt.contourf(DensityMax2D, cmap=plt.cm.Blues) #WORKS, but there are degenerated data + not enough data for interpolation. 
  #CS = plt.contourf(1E9*np.transpose(DensityMax2D[:,0]),1E9*DensityMax2D[:,1],np.log10(DensityMax2D), cmap=plt.cm.RdBu_r) #WORKS
  plt.colorbar(CS)
  #plt.contourf(wavelength_uniq, Efield1_uniq, DensityMax2D) #, label=r'$N_{exc}^{final}$') #, I='+str(libKeldysh.FieldToIntensity(Efield1_value))+r' W/m$^2$')
  del PlottingSinglePulses_Efield
      
  #ax1.legend(loc=3, borderaxespad=0., bbox_to_anchor=(1.05, 0.05))
  plt.xlabel(r'Laser wavelength (um)')
  plt.ylabel(r'Field amplitude (cell number)')
  #plt.ylabel(r'Final density (m$^{-3}$)')
  #plt.ylabel(r'Final density (m$^{-3}$)')
  #plt.xlim((100e-9,3100e-9))
  #ax1.set_xscale('log')
  #ax1.set_yscale('log')

  # Saving the plot
  plt.tight_layout()
  plt.savefig('Metaplot-SinglePulse-SingleColor-EffectOfWavelengthAndField.eps')
  plt.savefig('Metaplot-SinglePulse-SingleColor-EffectOfWavelengthAndField.png')
  if(ShowFinalPlot==1):
    plt.show()

  del PlottingSinglePulses, DensityMax2D
#}}}

