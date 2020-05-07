#!/usr/bin/env python
#-*- coding: utf-8 -*-

# Copyright (C) 2013-2020 T. J.-Y. Derrien
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

## @package libAtomicUnits
# Manage conversion between atomic units (AU) and International System (SI)

from numpy import vectorize
from scipy.constants import c, epsilon_0, mu_0, pi, e, m_e, h, hbar

BohrRadius=5.2917721067e-11 #m
Hartree_Energy = 27.21138602 #eV
c_AU = 137.035999139 
## Converts a charge density from AU to SI (m^-3).
def Density_AU_to_SI(density_AU):#{{{
  return density_AU * BohrRadius**(-3)
#}}}

## Converts a volume from AU to SI (m^{-3}).
def Volume_AU_to_SI(volume_AU):#{{{
  return volume_AU * BohrRadius**3
#}}}

## Converts time from AU to SI (seconds)
def Time_AU_to_SI(time_AU):#{{{
  return time_AU * hbar / e / Hartree_Energy
#}}}

## Converts time from SI to AU (?)
def Time_SI_to_AU(time_SI):#{{{
  return time_SI / (hbar / e / Hartree_Energy)
#}}}

## Converts length from SI (meters) to AU (Bohr radius)
def length_SI_to_AU(length_SI): #{{{
  return length_SI/BohrRadius
#}}}

## Converts length from AU (Bohr radius) to SI (meters) 
def length_AU_to_SI(length_AU): #{{{
  return length_AU*BohrRadius
#}}}

## Converts Electric field (Hartree / Bohr) to SI (V/m).
def Field_AU_to_SI(Field_SI):#{{{
  return Hartree_Energy / BohrRadius * Field_SI
#}}}

## Converts Electric field SI (V/m) to (Hartree / Bohr).
def Field_SI_to_AU(Field_AU):#{{{
  return (Hartree_Energy / BohrRadius) ** (-1.) * Field_AU
#}}}

def Energy_Hartree_to_eV(Energy_Ha):#{{{
  return Energy_Ha * Hartree_Energy
#}}}
def Energy_eV_to_Hartree(Energy_eV):#{{{
  return Energy_eV / Hartree_Energy
#}}}

## Converts electrical current from AU to SI
def Current_AU_to_SI(Current_AU):
    return Hartree_Energy/BohrRadius/e * Current_AU

## Converts electrical current from SI to AU
def Current_SI_to_AU(Current_SI):
    return Current_SI * e * BohrRadius / Hartree_Energy

Current_AU_to_SI = vectorize(Current_AU_to_SI)
Current_SI_to_AU = vectorize(Current_SI_to_AU)
Energy_Hartree_to_eV = vectorize(Energy_Hartree_to_eV)
Energy_eV_to_Hartree = vectorize(Energy_eV_to_Hartree)
Density_AU_to_SI = vectorize(Density_AU_to_SI)
Volume_AU_to_SI  = vectorize(Volume_AU_to_SI)
Time_AU_to_SI    = vectorize(Time_AU_to_SI)
length_SI_to_AU  = vectorize(length_SI_to_AU)
