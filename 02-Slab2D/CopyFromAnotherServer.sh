#!/bin/bash
set -e

if [ -z $1 ]; then
	echo "Usage: Localize yourself on the TARGET server. "
	echo "And run ./CopyCluster2Cluster.sh <draco2eos | eos2draco>"
	exit
fi

#DryRun="--dry-run"
DryRun=""

DracoPrev=/ptmp
EOSPrev=/scratch

CommonPath=tjyderri/octopus/02-Slab2D/Si-4x4x6-k-1x1x1-L24-Mesh0.434637008855-Mesh0.340150702582/lambda800e-9-tau10e-15-field9E9-CEP0.0/lambda400e-9-tau10e-15-field9E9-CEP0.0-delay0.0/ #Dont forget the slash. 

#tjyderri/octopus/02-Slab2D/Si-4x4x6-k-1x1x1-L12-Save/

case $1 in 
	draco2eos )
		Target=${EOSPrev}/${CommonPath}
		Source=${DracoPrev}/${CommonPath}
		SourceServerAddress=tjyderri@draco.mpcdf.mpg.de
	;;
	eos2draco )
		Source=${EOSPrev}/${CommonPath}
		Target=${DracoPrev}/${CommonPath}
		SourceServerAddress=tjyderri@eos01.bc.rzg.mpg.de
	;;
	*)
	echo "This entry is not taken into account."
	;;
esac

echo "Source Directory: ${Source}"
echo "Target Directory: ${Target}"

mkdir -p ${Target}

rsync ${DryRun} -avhre ssh --progress ${SourceServerAddress}:${Source} ${Target}

echo "Data were copied into: ${Target}"
