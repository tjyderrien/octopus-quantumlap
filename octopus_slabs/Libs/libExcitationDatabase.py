#!/usr/bin/env python3
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2021 T. J.-Y. Derrien
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

## @package libExcitationDatabase
# Library of functions for collecting the generated octopus simulations from subfolders
# Builds a table with laser parameters, final densities, and final energies after simulation. 

# import libDatabase
# from importlib import import_module
# libDatabase = import_module("")

from   octopus_slabs.Libs.libLogging         import init_logger
import numpy as np
import spp_extended_theory.Libs.libDatabase as libDatabase
from   octopus_slabs.Libs.libPlotting import ListAllOutput

Header="[libExcitationDatabase] "
logger = init_logger(__name__, verbose=True) #"plotFinalQuantities")

## Calls for the construction of excitation database and manages exceptions. 
def GetExcitationDatabase(path='*/*/*/td.general/', DatabaseName = 'Database.csv'):
  #DatabaseName = 'Database.csv'
  ## DEBUG PURPOSE ONLY
  # path="Si-1x1x1-k-24x24x1-Mesh0.434637008855-Mesh0.434637008855-T0e0/lambda800e-9-tau30E-15-field5E8-CEP0.0-StH2/lambda1600e-9-tau30E-15-field9E9-CEP0.0-delay0e-15-StH/td.general/"
  FromScratch = False #False: we try to re-use the existing DB. #True: we rebuild the database. 
  fileheader = "#T. J.-Y. Derrien, HiLASE Centre, derrien@fzu.cz.\n#Database v2.2.1\n#JobID\tRunDate\tkgrid\trealmesh\tlambda1\ttau1\tEfield1\tPolarizationAngle\tCEP1\tlambda2\ttau2\tEfield2\tCEP2\tdelay\tCrossPolar\tTheoryLevel\tInducedFields\tProbeFieldAmplitude_AU\tProbeFieldDelay_s\tFunctional\tFrozenHxc\tMoveIon\tNeMax\tNeFinal\tVirtualOscillations\tFinalCurrent\tMaxCurrent_AU\tMaxJx\tMaxJy\tMaxJz\tEnergy_eV\tEigenvalues_eV\tHartree_Energy_eV\tExchangeEnergy_eV\tEnergy_afterpulse_eV\tEigenvalues_afterpulse_eV\tHartree_Energy_afterpulse_eV\tExchangeEnergy_afterpulse_eV\tPoyntingEnergyFinal_eV\tEext_peak_V_m\tN_exc_afterpulse\tWpi\tWpi_afterpulse\tSTH\tRiseCycles\tFinalTime\tCaptureTime"

  logger.info(Header+"Attempting to load previous results...")

  try:
    logger.info("** Info: just loading existing database. ")
    Database = np.loadtxt(DatabaseName, skiprows=1)
    logger.info(Database)
    logger.info("Previous database was found. Attempting to complete it.")
    FromScratch=False #No need to rebuild database from scratch
#    ListAllOutput(path, DatabaseName, FromScratch)
  except:
    logger.warning(Header+"Previous database could NOT be found. Rebuilding the database from scratch.")
    FromScratch = True #We need to rebuild database from scracth

  if(FromScratch): #rebuilds database from scratch
    FreshData, NewDatabase = ListAllOutput(path, DatabaseName, FromScratch)
    logger.info(Header+"Save the table to a new file...")
    libDatabase.ExportToTxt(NewDatabase, DatabaseName, fileheader)
  else: 
    logger.warning("Entered case with FromScratch=False. Newly found entries may not be saved. ")
    # ListAllOutput(path, DatabaseName, FromScratch)
    NewDatabase=Database
  return NewDatabase
