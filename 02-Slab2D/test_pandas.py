#!/usr/bin/env python2
#-*- coding: utf-8 -*-

# Copyright (C) 2016-2020 T. J.-Y. Derrien
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

## @package test_pandas
# It seems one can manage data much more easily using pandas. 
# Lets if this saves time.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Loading database
data = pd.read_csv("Database.csv", sep="\t", header=2) #pandas can use heterogeneous data, and manages header automatically (then program of fdominec is over)

# Selecting subsets for repeating PRL Fig. 4. 
## Selecting wavelengths and pulse duration
data = data[ 
    ( data["lambda1"]==3200e-9 ) & 
    ( data["tau1"]==20e-15 ) & 
    ( data["CEP1"]==0 ) &
    ( data["STH" ]==3 )
    ] 

print(data)

subdata=data[["lambda1", "Efield1", "NeFinal"]]

# Plotting these data
subdata.plot()
#data.plot_scatter(x="wavelength1", y="field1", alpha=0.5)
plt.show()
