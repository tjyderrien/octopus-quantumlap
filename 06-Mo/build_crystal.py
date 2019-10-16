#!/bin/python

from ase.build import bulk
interatomic_dist = 3.147 #Angstroms #https://materials.springer.com/isp/crystallographic/docs/sd_0554360
Mo=bulk("Mo", "bcc", interatomic_dist, cubic=True)
print Mo.cell
