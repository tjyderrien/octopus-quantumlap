#!/bin/bash
## WARNING: it seems ${BASH_SOURCE} generates serious mistakes, such as bugs like SCP file transfers get blocked.

SourcePath="${BASH_SOURCE[0]}"
LocalFolder="$( cd "$(dirname "${SourcePath}" )" && pwd)" #Should work on any computer but can generate many problems, such as stopping SCP transfers.
source ${LocalFolder}/../machines.env 
CheckComputerName
LoadEnvLibraries

case "${ComputerName}" in 
	${LenovoTibo})
		export QuantumLaPsources="/home/hilase/Documents/octopus-slabs"
		export QuantumLaPruns="/home/hilase/Documents/octopus"
	;;

	${PrometheusLogin}*)
		export QuantumLaPruns=/net/scratch/projects/prace/pr2c0043
		export QuantumLaPsources=$HOME/octopus-slabs
	;;

	${SalomonLogin}*)
		export QuantumLaPsources=/home/tjyderrien/octopus-slabs
		export QuantumLaPruns=/scratch/work/project/open-16-27/octopus-slabs.data
	;;

        ${BIATRI}*)
		export QuantumLaPsources="/home/thibault.derrien/Documents/octopus-slabs"
		export QuantumLaPruns="/home/thibault.derrien/Documents/octopus.data"
	;;
        *)
            ## DO NOTHING
        ;;


#		export QuantumLaPsources="~/Documents/octopus-slabs"
#        export QuantumLaPruns="~/Documents/octopus.data"
esac

