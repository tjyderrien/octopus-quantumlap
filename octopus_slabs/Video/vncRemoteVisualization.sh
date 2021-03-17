#!/bin/bash
# echo "Usage: ./script.py <server|client> <Remote Port Number>"
Commuter=$1
GraphPortNumber=$2

if [ -z $1 ]; then
	echo "Usage: ./script.py <server|client> <Remote Port Number> [<Remote machine name>]"
	exit
fi

source ../02-Slab2D/machines.env
CheckComputerName

#Defining clients and servers
echo "Computer name: $ComputerName"
# if [ ${ComputerName} == ${EOS}* ]; then
#Commuter=client
#fi
#if [ ${ComputerName} == ${AsusThibs} ]; then
#	Commuter=client
#	echo "Client."
#fi

#GraphPortNumber=$1
# ServerAddress=eos01.bc.rzg.mpg.de #included in machines.env
# LoginName=tjyderri #also managed by machines.env
# GateWayAddress=gate.rzg.mpg.de #managed by machines.env

if [ ${Commuter} == "client" ]; then
	echo "Local computer is client."
	if [ -z $3 ]; then
		echo "Please indicate name of the machine on 3rd argument."
		echo "Possibilities: $Drako, $Hydra, $EOS"
		exit -1
	else
        	echo "Remote server is: $3"
		vncviewer -via ${LoginName}@${GateWayAddress} ${ServerAddress}:${GraphPortNumber}
	fi
fi

if [ ${Commuter} == "server" ]; then
	echo "Local computer is server."
	echo "Loading modules..."
#	source loadModules.sh
#	ModuleForEOS
	echo "Start the VNC server."
	vncserver
fi

# echo "This script is just a documentation."
# echo "On the REMOTE machine:"
# echo "Load the modules given in VisualizationWorkshop/loadModules.sh"
# echo "vncserver"
# echo "New 'X' desktop is eos01:1"

# echo "On the LOCAL machine:"
# echo "sudo apt-get install vncviewer"
# echo "vncviewer -via tjyderri@gate.rzg.mpg.de eos01.bc.rzg.mpg.de:1"
# echo "vncviewer -via tjyderri@gate.rzg.mpg.de draco.mpcdf.mpg.de:1
