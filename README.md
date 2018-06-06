# Description

Private repository made to develop Octopus simulations for QuantumLaP project. 

WARNING: the repository spp-extended-theory is necessary to run the routines. Download the program here: git clone git@bitbucket.org:tjyderrien/spp-extended-theory.git.  

WARNING2: before having a successful usage, you might have to link some libraries by hand with symbolic links. #TODO: automake should do that for us someday. Or simply linking libraries properly in Python. 

# Documentation

dOxygen documentation is available in ./doc. Generate using generate.sh. 

# Languages

The routines are a combination of Python and Bash. 
Of course, it may not be completely optimum. 

# List of folders

* Papers/: contains the future possible papers

* FHI-aims/: Molecular dynamics software to build up the lattice structures: FHI-aims. 

* The other folders (00..xx/ contain octopus simulation files).

* 00-SelfMadeFiles: Files made from the octopus tutorials. 

* 01-Example given obtained from Umbierto

* 02-Slab2D: examples of graphene, Si-H samples

* 03-Bulk: everything you need to compute bulks. 

* Samples: preparation of bulk and slabs samples. 
