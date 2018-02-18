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

import libDatabase #have to separate to avoid conflicts
from libPlotting import ListAllOutput

Header="[libExcitationDatabase] "

def GetExcitationDatabase(path):
#  path='*/*/td.general/'
  DatabaseName = 'Database.csv'

  print Header+"Attempting to load previous results..."
  try:
    PlottingTable = loadtxt(DatabaseName, delimiter="\t")
  except:
    print Header+"Not found. Rebuilding the index."
    PlottingTable = ListAllOutput(path)

    print Header+"Save the table to a file..."
    try:
      libDatabase.ExportToTxt(PlottingTable, DatabaseName)
    except:
      print "** Error when writing the database.csv."
      exit()

  return PlottingTable
