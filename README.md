# Code author: 
Thibault J.-Y. Derrien <derrien@fzu.cz>
HiLASE Centre, Institute of Physics (Prague, Czech Academy of Sciences)
Max Planck Institute for Structure and Dynamics of Matter (Hamburg, Germany)

# Description

Private repository made to develop Octopus simulations for: 
- MSCA-IF-2015/QuantumLaP (prof. Bulgakova / T. Derrien)
- CZ/BIATRI (prof. Bulgakova)
- MSCA-RISE-2018/ATLANTIC (T. Derrien)

WARNING: the repository spp-extended-theory is necessary to run the routines. Download the program here: git clone git@bitbucket.org:tjyderrien/spp-extended-theory.git.  

WARNING2: before having a successful usage, you might have to link some libraries by hand with symbolic links. #TODO: automake should do that for us someday. Or simply linking libraries properly in Python. 

# Versioning information
Versioning strategy is as follows. 
1. Each new format/structure of database is a major change, and therefore defines a new version v0.1, v0.2, etc.
2. New feature or retro-active modification of the code is a minor change, and defines a sub-version v0.1.2, etc.

Since this logic would induce already be at v13 when enabling version tagging, it was decided to adopt this less hurry scheme.

# Databases

Several databases are available. For now, they are separated into :
* Material, i.e., periodic lattice definition, 
* Ground state k-space grid,
* Ground state real-space grid,
* Groud state Fermi temperature,
* Employed xc functional.

# Documentation
dOxygen documentation is available in ./doc. Generate using doc/generate.sh. 

# Direction/Languages

The routines are a combination of Python and Bash. 
Of course, it may not be completely optimum and transition to ASE (Python) would be desirable in future.

# List of folders
Warning: avoid putting simulations into the git tree. Since generated data can be huge (> TB), 

* FHI-aims/: Some examples for FHI-aims are provided. Molecular dynamics software to build up the lattice structures: FHI-aims. This folder shall be removed at some point. 

* 00-SelfMadeFiles: Files made from the octopus tutorials. 

* 01-Examples given by Umbierto de Giovannini. 

* 02-Slab2D: examples of graphene, Si-H samples

* 03-Bulk: everything you need to compute Si bulks. 

* xx-Material: many materials start to be available. Routines to compute them in series are given in these folders.

* Samples: preparation of bulk and slabs samples. 

* Pulses: several types of laser pulses have been introduced into the routines: sin**2, sin**4, soft top-hat (STH), Gaussian (disadvised).

* doc: in this folder, automatic dOxygen documentation can be generated easily.
