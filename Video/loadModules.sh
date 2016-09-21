#!/bin/bash
# TODO: make this machine free. 
# source ../02-Slab2D/machines.env

ModuleForEOS () {
	module purge
	module load remotevis
	module load ffmpeg/2.7
	module load python27/python/2.7
<<<<<<< HEAD
	module load python27/scipy/2015.01
=======
>>>>>>> 37c8b3b9c5e25f5bb5b61d6e953b0a986a99de78
	module load visit/2.9
}

ModuleForDraco() {
	module purge
	module load ffmpeg/3.0
	module load anaconda/2
	module load visit/2.10 
}



echo "For making a video using the Py script: visit --cli --nowin"
