#!/bin/python

from ase.build import bulk
interatomic_dist = 3.147 #Angstroms #https://materials.springer.com/isp/crystallographic/docs/sd_0554360
SiO2=bulk("SiO2") #cubic=True)
print(SiO2.cell)
