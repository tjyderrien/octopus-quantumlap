#!/bin/bash
LocalFolder=$(dirname $BASH_SOURCE) #Should work on any computer.
source $LocalFolder/../machines.env
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

	*)
		export QuantumLaPsources="~/Documents/octopus-slabs"
                export QuantumLaPruns="~/Documents/octopus.data"
	;;
esac



