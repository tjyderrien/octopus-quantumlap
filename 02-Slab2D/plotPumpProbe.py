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

# @module plotPumpProbe
# Provides a set of functions aimed at analyzing the optical properties of solids from the total_current. 

import numpy as np
from libAtomicUnits import *
from libPlotting import *

# 1. Select 2 runs, 

path_run1 = "lambda800e-9-tau20E-15-field1E10-CEP0.0-SinSquared2-Probe0.1-ProbeDelay1e-15/lambda400e-9-tau20E-15-field0E0-CEP0.0-delay0.5e-15-SinSquared-InducedFields/td.general/total_current" #PUMP 
path_run2 = "lambda800e-9-tau20E-15-field1E10-CEP0.0-SinSquared2-Probe1.0-ProbeDelay1e-15/lambda400e-9-tau20E-15-field0E0-CEP0.0-delay0.5e-15-SinSquared-InducedFields/td.general/total_current" #PUMP+PROBE RUN

# 2. Capture their total_current. 
current1=np.loadtxt(path_run1, dtype='str', delimiter="\t", skiprows=4)
current2=np.loadtxt(path_run2, dtype='str', delimiter="\t", skiprows=4)
# Considered using ASE here, but it is poorly developed for TDDFT. 

# 3. Identify the common temporal intervals. 
current1_time_AU, current1_X_AU, current1_Y_AU, current1_Z_AU = GetTotalCurrent(current1)
current2_time_AU, current2_X_AU, current2_Y_AU, current2_Z_AU = GetTotalCurrent(current2)

current12_time_AU_min = np.max(np.min(current1_time_AU), np.min(current2_time_AU))
current12_time_AU_max = np.min(np.max(current1_time_AU), np.max(current2_time_AU))

# 4. Interpolate the signals if necessary (seems not required here). 
# 5. Prepare Yabana's case: J(pump+probe) - J_pump = J_probe. 
current_probe_X_AU=current2_X_AU-current1_X_AU
current_probe_Y_AU=current2_Y_AU-current1_Y_AU
current_probe_Z_AU=current2_Z_AU-current1_Z_AU

# 6. Fourier Transform the J_probe. 
# 7. Extract the conductivity out of it. 
# 8. Deduce the dielectric property. 
