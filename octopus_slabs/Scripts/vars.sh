#!/bin/bash
## WARNING: it seems ${BASH_SOURCE} generates serious mistakes, such as bugs like SCP file transfers get blocked.

# LocalFolder=$(dirname ${BASH_SOURCE} ) #Should work on any computer but can generate many problems
# source ${LocalFolder}/../machines.env 

SourcePath="${BASH_SOURCE[0]}"
# LocalFolder="$( cd "$(dirname "${SourcePath}" )" && pwd)" #Should work on any computer but can generate many problems, such as stopping SCP transfers, or interfering with VNC and X11 servers.
# echo "Debug: \$SourcePath: $SourcePath"
# LocalFolder=$( dirname $SourcePath )

# Safer version given by ChatGPT
# LocalFolder="$(cd "$(dirname -- "$SourcePath")" >/dev/null 2>&1 && pwd -P)"

# Safer version given by ChatGPT with providing context as source vars.sh

pushd "$(dirname -- "$SourcePath")" >/dev/null 2>&1
LocalFolder="$(pwd -P)"
popd >/dev/null 2>&1

# echo "Debug: \$LocalFolder: $LocalFolder"
source ${LocalFolder}/../machines.env 
# source ../machines.env

# Because of SSH transfers, this file must not emit any return. 

CheckComputerName > /dev/null
echo "[octopus_slabs: vars.sh] Computer is : $ComputerName"
LoadEnvLibraries  > /dev/null

case "${ComputerName}" in 
	${LenovoTibo})
		export QuantumLaPsources="/home/hilase/Documents/octopus-slabs"
		export QuantumLaPruns="/home/hilase/Documents/octopus"
		echo "Stored variables for local laptop."
	;;

	*${PrometheusLogin}*)
		export QuantumLaPruns="/net/scratch/projects/prace/pr2c0043"
		export QuantumLaPsources="$HOME/octopus_slabs"
	;;

	${SalomonLogin}*)
		export QuantumLaPsources=/home/tjyderrien/octopus-slabs
		export QuantumLaPruns=/scratch/work/project/open-20-5/octopus-slabs.data
	;;

        ${BIATRI})
		export QuantumLaPsources="/home/thibault.derrien/Documents/octopus-slabs"
		export QuantumLaPruns="/home/thibault.derrien/Documents/octopus.data"
	;;

${biatri2})
	export QuantumLaPsources="/home/tderrien/Documents/octopus-slabs"
	export QuantumLaPruns="/home/tderrien/Documents/octopus.data"
	echo "[octopus-slabs] Prepared QuantumLaP variables for BIATRI2."
	;;

	${kastanienallee})
		export QuantumLaPsources="/home/thibault/Documents/Codes/octopus-slabs"
		export QuantumLaPruns="/home/thibault/Documents/Codes/octopus.data"

	;;

        ${IBM})
                export QuantumLaPsources="/home/tjyderrien/octopus-slabs"
                export QuantumLaPruns="/home/tjyderrien/octopus.data"

        ;;

        ${BarboraLogin} | ${BarboraRuns})
                export QuantumLaPsources="/home/tjyderrien/octopus-slabs"
                export QuantumLaPruns="/scratch/project/open-20-5/octopus-slabs.data"
        ;;

        ${navigator_login} | ${navigator}* )
                echo "[vars.sh: ] Definition QuantumLaP variables for Navigator."
                export QuantumLaPsources="/veracruz/home/p/pr2c0043/octopus-slabs/"
                export QuantumLaPruns="/tormenta/f/frecuencia/octopus-slabs.data/"
        ;;

${NvidiaDGX})
	 echo "[vars.sh: ] Definition QuantumLaP variables for Masaryk Uni Nvidia DGX."
	 export QuantumLaPsources=/storage/brno12-cerit/home/derrit01/octopus-slabs/
         export QuantumLaPruns=/storage/brno12-cerit/home/derrit01/octopus.data/

	;;
        

        *)
            ## DO NOTHING
	    echo "[octopus_slabs: vars.sh] Case not covered. Computer: ${ComputerName}"
        ;;


#		export QuantumLaPsources="~/Documents/octopus-slabs"
#        export QuantumLaPruns="~/Documents/octopus.data"
esac

source ${QuantumLaPsources}/Scripts/functions.sh

