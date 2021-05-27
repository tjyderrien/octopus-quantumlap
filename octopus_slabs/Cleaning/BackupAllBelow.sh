#!/bin/bash
set -e
source $QuantumLaPsources/octopus_slabs/machines.env
source $QuantumLaPsources/Scripts/functions.sh
CheckComputerName
LoadEnvLibraries

source $QuantumLaPsources/Cleaning/TarResults.sh

backup_StH_fixed > backup-StH.log
backup_SinSquared > backup-SinSquared.log
backup_SinFourth > backup-SinFourth.log
