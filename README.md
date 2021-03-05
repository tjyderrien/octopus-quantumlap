# Code author: 
Thibault J.-Y. Derrien <derrien@fzu.cz>
HiLASE Centre, Institute of Physics (Prague, Czech Academy of Sciences)
Max Planck Institute for Structure and Dynamics of Matter (Hamburg, Germany)

# Coding standards

To converge and merge with Atomic Environment Simulation (ASE) in future, we aim at adopting coding standards detailed 
here: <https://wiki.fysik.dtu.dk/ase/development/python_codingstandard.html>

# Description

Private repository made to develop Octopus simulations for: 
- MSCA-IF-2015/QuantumLaP (prof. N. M. Bulgakova / T. Derrien)
- CZ/BIATRI (prof. Bulgakova)
- MSCA-RISE-2018/ATLANTIC (T. Derrien)

WARNING: the repository "spp_extended_theory" is necessary to run the routines. Download the program here: 
git clone git@bitbucket.org:tjyderrien/spp_extended_theory.git.  

WARNING2: if you meet necessity of employing symbolic linking, write to developer! Solution for this is now understood.

# Versioning information
Versioning strategy is as follows. 
1. Each new format/structure of database is a major change, and therefore defines a new version v0.1, v0.2, etc.
2. New feature or retro-active modification of the code is a minor change, and defines a sub-version v0.1.2, etc.

# Databases

Several databases generated using this code are available. 
Since v2.0, there is only 1 database per material. Each database contains all the TD-runs made with various functionals,
various k-grid precision, mesh precision, pulses shapes, etc.  

# Documentation
Automatic dOxygen documentation can be generated in ./doc using doc/generate.sh. "Dot" picture utility is necessary to 
generate the automatic diagrams.  

Tutorials, documentation, books, etc can be found on http://www.quantumlap.eu/

# Direction/Languages

The routines are a combination of Python and Bash.
Bash is more long-lasting investment in time. 
Transition to ASE (Python) would be desirable in future.
Python is more trendy and more modular, and helps to prepare transition to ASE.  

# List of folders
Warning: don't upload simulations results into the git tree. Generated data can be huge (> TB) and Git is NOT meant 
to contain binary data. One thing exists from Git to host binary data: Git-LFS.

At the moment, folders starting with a number contain materials & files to start heavy DFT and TDDFT computations. 

Folders below octopus_slabs/octopus_slabs/ contain libraries and post-treatment scripts for analysis.   

* 00-SelfMadeFiles: Files made from the octopus tutorials. 

* 01-Examples given by Umbierto de Giovannini. 

* 02-Slab2D: examples of graphene, Si-H samples

* 03-Bulk: everything you need to compute Si bulks. 

* xx-Material: many materials start to be available. Routines to compute them in series are given in these folders.

* Samples: preparation of bulk and slabs samples. 

* Pulses: several types of laser pulses have been introduced into the routines: sin**2, sin**4, soft top-hat (STH), Gaussian (disadvised).

* doc: in this folder, automatic dOxygen documentation can be generated easily.
