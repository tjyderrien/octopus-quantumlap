#!/bin/bash
set -e
source $QuantumLaPsources/machines.env
source $QuantumLaPsources/02-Slab2D/functions.sh
CheckComputerName
LoadEnvLibraries

source $QuantumLaPsources/Cleaning/TarResults.sh

backup_StH_fixed > backup-StH.log
backup_SinSquared > backup-SinSquared.log
backup_SinFourth > backup-SinFourth.log
