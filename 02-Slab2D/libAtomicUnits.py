#!/usr/bin/env python
#-*- coding: utf-8 -*-

## @module libAtomicUnits
## Manage conversion between atomic units (AU) and International System (SI)

BohrRadius=5.2917721067e-11 #m

## Converts a charge density from AU to SI (m^-3).
def Density_AU_to_SI(density_AU):#{{{
  return density_AU * BohrRadius**(-3)
#}}}

## Converts a volume from AU to SI (m^{-3}).
def Volume_AU_to_SI(volume_AU):#{{{
  return volume_AU * BohrRadius**3
#}}}