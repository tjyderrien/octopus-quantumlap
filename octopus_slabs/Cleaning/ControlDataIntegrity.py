#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2023-2023 T. J.-Y. Derrien
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

## @package ControlDataIntegrity
# Define shared functions for plotting routines and database construction

from octopus_slabs.Libs.libString import contains_binary_data

## Controls the integrity of LOG files in the amount of accumulated octopus data
def CheckLOGfileIntegrity(path="*/*/*/td.general/")
    # import time
    for FolderName in glob.glob(path):
         # start_time = time.time()
         LOG_newest = max(glob.glob(FolderName+'../octopus*.log'), key=os.path.getmtime)
         if(contains_binary_data(LOG_newest)):
             print(LOG_newest)
         # else:
         #    print(FolderName, LOG_newest)
         # end_time = time.time()
         # elapsed_time = end_time - start_time
         # print(f"Elapsed time: {elapsed_time} seconds")
    return 0
 




