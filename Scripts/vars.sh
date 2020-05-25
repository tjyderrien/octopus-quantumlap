#!/bin/bash

source ../machines.env
CheckComputerName

case "${ComputerName}" in 
	${LenovoTibo})
		export QuantumLaPsources="~/Documents/octopus-slabs"
		export QuantumLaPruns="~/Documents/octopus.data"
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

