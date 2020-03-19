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

# 1. Select 2 runs, 
# 2. Capture their total_current. 
# 3. Identify the common intervals. 
# 4. Interpolate the signals if necessary. 
# 5. Prepare Yabana's case: J(pump+probe) - J_pump = J_probe. 
# 6. Fourier Transform the J_probe. 
# 7. Extract the conductivity out of it. 
# 8. Deduce the dielectric property. 
