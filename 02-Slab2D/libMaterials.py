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


## @package libMaterials
# Puts there the properties of materials we are modeling

# 1. Select the material one want to model. This can be based on the name of the folder. 
# 2. Load the relevant properties for this material. 

def SiBulk(): #{{{
  Egap_SI = 2.56*e #[Waroquiers and Gonze, Phys. Rev. B 87, 075121 (2013)]
  meff = 0.2226    #[Laflamme Janssen and X. Gonze, Physical Review B 93, 205147 (2016)]
  OpticalIndex={
        400e-9: 5.568072776280323, 483e-9: 4.4028, 484e-9: 4.3964, 
        600e-9: 3.947144634525661,
        800e-9: 3.6924167231905, 1030e-9: 3.57807133869753, 1200e-9: 3.51933432518597, 1600e-9: 3.4826,
        2200e-9: 3.4548, 2600e-9: 3.4457, 3200e-9: 3.43438589865999}
  return Egap_SI, meff, OpticalIndex
