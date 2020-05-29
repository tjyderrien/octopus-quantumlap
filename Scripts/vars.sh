#!/bin/bash
## WARNING: don't call this program using "source", otherwise bugs like SCP file transfers get blocked.

LocalFolder=$(dirname $BASH_SOURCE) #Should work on any computer. 
$LocalFolder/../machines.env #BUG: Using "source" here kills transfers via scp or rsync 
CheckComputerName

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
		export QuantumLaPsources="~/Documents/octopus-slabs"
        export QuantumLaPruns="~/Documents/octopus.data"
	;;
esac

