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

# @module libExcitationDatabase
# Library of functions for collecting the generated octopus simulations from subfolders
# Builds a table with laser parameters, final densities, and final energies after simulation. 

import libDatabase 
from libPlotting import ListAllOutput
from numpy       import loadtxt

Header="[libExcitationDatabase] "

## Calls for the construction of excitation database and manages exceptions. 
def GetExcitationDatabase(path='*/*/td.general/', DatabaseName = 'Database.csv'):
#  path='*/*/td.general/'
  #DatabaseName = 'Database.csv'
  FromScratch = False #always try to re-use the existing DB. 
  print Header+"Attempting to load previous results..."
  try:
    print Header+"** Info: loading existing database. "
    Database = np.loadtxt(DatabaseName, skiprows=0)
    print Database
    print Header+"** Info: Previous database was found. Attempting to complete it."
    FromScratch=False
    libPlotting.ListAllOutput(path, DatabaseName, FromScratch)
  except:
    print Header+"** Warning: Previous database could NOT be found. Rebuilding the database from scratch."
    FromScratch = True
    FreshData, NewDatabase = ListAllOutput(path, DatabaseName, FromScratch)

    print Header+"Save the table to a new file..."
    libDatabase.ExportToTxt(NewDatabase, DatabaseName)
  return NewDatabase
