#!/usr/bin/env python
#-*- coding: utf-8 -*-

## @module libAtomicUnits
## Manage conversion between atomic units (AU) and International System (SI)

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
  return hbar / e / 27.21138505e0
#}}}

## Converts length from SI (meters) to AU (Bohr radius)
def length_SI_to_AU(length_SI): #{{{
  return length_SI/BohrRadius
#}}}

## Converts Electric field (Hartree / Bohr) to SI (V/m).
def Field_AU_to_SI(Field_SI):#{{{
  return Hartree_Energy / BohrRadius * Field_SI
#}}}
Density_AU_to_SI = vectorize(Density_AU_to_SI)
Volume_AU_to_SI  = vectorize(Volume_AU_to_SI)
Time_AU_to_SI    = vectorize(Time_AU_to_SI)
length_SI_to_AU  = vectorize(length_SI_to_AU)
