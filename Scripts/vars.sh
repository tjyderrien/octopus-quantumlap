#!/bin/bash
## WARNING: it seems ${BASH_SOURCE} generates serious mistakes, such as bugs like SCP file transfers get blocked.

# LocalFolder=$(dirname ${BASH_SOURCE} ) #Should work on any computer but can generate many problems
# source ${LocalFolder}/../machines.env 

SourcePath="${BASH_SOURCE}[0]"
LocalFolder="$( cd "$(dirname "${SourcePath}" )" && pwd)" #Should work on any computer but can generate many problems, such as stopping SCP transfers.
source ${LocalFolder}/../machines.env 

# Because of SSH transfers, this file must not emit any return. 

CheckComputerName > /dev/null
echo "[octopus-slabs: vars.sh] Computer is : $ComputerName"
LoadEnvLibraries  > /dev/null

case "${ComputerName}" in 
	${LenovoTibo})
		export QuantumLaPsources="/home/hilase/Documents/octopus-slabs"
		export QuantumLaPruns="/home/hilase/Documents/octopus"
		echo "Stored variables for local laptop."
	;;

	*${PrometheusLogin}*)
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
        ${BarboraLogin} | ${BarboraRuns})
                export QuantumLaPsources=/home/tjyderrien/octopus-slabs
                export QuantumLaPruns=/scratch/project/open-20-5/octopus-slabs.data
        ;;
        
        *)
            ## DO NOTHING
	    echo "[octopus-slabs: vars.sh] Case not covered. Computer: ${ComputerName}"
        ;;


#		export QuantumLaPsources="~/Documents/octopus-slabs"
#        export QuantumLaPruns="~/Documents/octopus.data"
esac

