#!/bin/bash

set -e

program=gnuplot

$program -e "field='1E9'" $QuantumLaPsources/02-Slab2D/plotTimeDensity-MultiParameters.gnu
$program -e "field='5E9'" $QuantumLaPsources/02-Slab2D/plotTimeDensity-MultiParameters.gnu
$program -e "field='1E10'" $QuantumLaPsources/02-Slab2D/plotTimeDensity-MultiParameters.gnu

$program -e "field='1E9'" $QuantumLaPsources/02-Slab2D/plotTimeEnergy-MultiParameters.gnu
$program -e "field='5E9'" $QuantumLaPsources/02-Slab2D/plotTimeEnergy-MultiParameters.gnu
$program -e "field='1E10'" $QuantumLaPsources/02-Slab2D/plotTimeEnergy-MultiParameters.gnu



