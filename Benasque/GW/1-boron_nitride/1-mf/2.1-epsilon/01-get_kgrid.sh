#!/bin/bash -l

module load berkeleygw/1.2
kgrid.x ./kgrid.inp ./kgrid.out ./kgrid.log &> out
