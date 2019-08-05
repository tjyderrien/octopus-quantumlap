#!/bin/python

from ase.build import surface
from ase.visualize import view
interatomic_dist = 3.147 #Angstroms #https://materials.springer.com/isp/crystallographic/docs/sd_0554360
#Mo=bulk("Mo", "bcc", cubic=True)
#print Mo.cell

Al=surface('Al', (1,1,1), 2)
Al.center(vacuum=5,axis=1)

view(Al)
