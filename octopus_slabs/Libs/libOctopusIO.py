#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""@package docstring
Module with functions to parse Octopus output files and load them into NumPy arrays and/or dictionaries.
Parsing of input files has not yet been implemented.

Module written by Andrés I. Bertoni (andresibertoni@gmail.com).
"""

import os
import copy
import pickle
import numpy as np

def save_pickle(data, output_path):
    """! Stores 'data' as a pickled binary file at 'output_path'.
    @param data Object to be stored as binary (e.g., python data structure or NumPy array).
    @param output_path Path of the file to be saved (output binary). It is recommended to include a .pkl extension.
    @return Does not return anything. But a binary file is writen instead.
    """
    picklefile = open(output_path,"wb")
    pickle.dump(data, picklefile)
    picklefile.close()

def load_pickle(pkl_path):
    """! Loads pickled object from binary file at 'pkl_path'.
    @param pkl_path Path of the pickled object file to be loaded (binary file with -probably- .pkl extension).
    @return Returns the loaded object (e.g., python data structure or NumPy array).
    """
    picklefile = open(pkl_path,"rb")
    data = pickle.load(picklefile)
    picklefile.close()
    return data

def load_k_paths(path_to_kfiles, return_uncat=False, u_to_idx={"_x":0, "_y":1, "_z":2}):
    """! Recognizes the paths of the k-points output files located in 'path_to_k_files'.
    @param path_to_kfiles Path to the folder with all the k-points output files.
    @param return_uncat If the uncategorized paths should be returned or not. Default value is False.
    @param u_to_idx Dictionary indicating which filename-ending corresponds to which index of spatial direction. Default value is {"_x":0, "_y":1, "_z":2}; those are the three Cartesian spatial directions.
    @return k_indices Array with the indices of the k-points that were successfully loaded in all spatial directions.
    @return k_paths A matrix containing, for each one of the spatial directions, the paths to the k-points files, sorted according to the k_index array.
    @return uncat_paths (return_uncat=True) List of paths to files in the same folder ('path_to_kfiles') that do not correspond to k-points output files or simply could not be categorized.
    """
    k_paths = [[] for idx in range(len(u_to_idx.keys()))]
    k_indices = copy.deepcopy(k_paths)
    uncat_paths = []
    # (separates paths based on direction index)
    for k_path in os.listdir(path_to_kfiles):
        u = k_path[-2:]
        if u in u_to_idx.keys():
            idx = u_to_idx[u]
            k_paths[idx].append(k_path)
            # expected path format: 'label.k(index)_u'
            k_index = int(k_path.replace(u,"").split(".")[-1][1:])
            k_indices[idx].append(k_index)
        else: uncat_paths.append(k_path)
    # (selects k indices available in three dimensions)
    ubiqui_k_indices = list(set(k_indices[0]) & set(k_indices[1]) & set(k_indices[2]))
    ubiqui_k_indices = np.array(ubiqui_k_indices) ; ubiqui_k_indices.sort()
    for u in u_to_idx.keys():
        idx = u_to_idx[u] ; u_str = u.replace("_","")
        ubiqui_mask = np.isin(k_indices[idx], ubiqui_k_indices)
        # Improve:
        #if not np.all(ubiqui_mask): print(f"Not all k-points are available for the {u_str} direction!")
        k_paths[idx] = np.array(k_paths[idx])[ubiqui_mask]
        k_indices[idx] = np.array(k_indices[idx])[ubiqui_mask]
    # (sorts paths based on k indices)
    for idx in u_to_idx.values():
        order = np.argsort(k_indices[idx])
        k_paths[idx] = np.array(k_paths[idx])[order].tolist()
    k_paths = np.array(k_paths)
    k_indices = ubiqui_k_indices
    # (returns output)
    if not return_uncat: return k_indices, k_paths
    else: return k_indices, k_paths, uncat_paths

def extract_length_dipole_elements(path_to_kfiles, u_to_idx={"_x":0, "_y":1, "_z":2}):
    """! Loads the dipole matrix elements (in length gauge) from the k-points output files located in 'path_to_k_files' into a NumPy array.
    @param path_to_kfiles Path to the folder with all the k-points output files.
    @param u_to_idx Dictionary indicating which filename-ending corresponds to which index of spatial direction. Default value is {"_x":0, "_y":1, "_z":2}; those are the three Cartesian spatial directions.
    @return k_indices Array with the indices of the k-points that were successfully loaded in all spatial directions.
    @return length_dipole_elements NumPy array with the loaded dipole matrix elements (in length gauge); it has shape (n_kpoints, n_dims, n_states, n_states), where n_kpoints is the number of available/recognized k-points, n_dims is the number of spatial dimensions in 'u_to_idx' and n_states is the number of 1-electron states available in the output.
    """
    k_indices, k_paths = load_k_paths(path_to_kfiles, u_to_idx=u_to_idx)
    n_dims, n_kpoints = k_paths.shape
    n_states = np.loadtxt(os.path.join(path_to_kfiles, k_paths.flatten()[0]), skiprows=3).shape[0]
    length_dipole_elements = np.zeros((n_kpoints, n_dims, n_states, n_states))
    for kdx in range(n_kpoints):
        for udx in range(n_dims):
            k_path = os.path.join(path_to_kfiles, k_paths[udx][kdx])
            length_dipole_elements[kdx,udx,...] = np.loadtxt(k_path, skiprows=3)
    return k_indices, length_dipole_elements

def store_length_dipole_elements(path_to_kfiles, output_path=".", kindx_pkl_fname="k_indices.pkl", lendipel_pkl_fname="length_dipole_elements.pkl", u_to_idx={"_x":0, "_y":1, "_z":2}):
    """! Extracts data from the k-points output files located in 'path_to_k_files' and store it as two separate pickled (binary) files, each containing a NumPy array: (1) with the extracted k-points indices and (2) with the dipole matrix elements (in length gauge).
    @param path_to_kfiles Path to the folder with all the k-points output files.
    @param output_path Path to the folder where the output pickled (binary) files should be stored. The default value is the current directory of the function execution (".").
    @param kindx_pkl_fname Filename for the pickled file with the extracted k-points indices (1). Default value is 'k_indices.pkl'.
    @param lendipel_pkl_fname Filename for the pickled file with the dipole matrix elements (in length gauge) (2). See the 'extract_length_dipole_elements' function in the same module for extended documentation. Default value is 'length_dipole_elements.pkl'.
    @param u_to_idx Dictionary indicating which filename-ending corresponds to which index of spatial direction. Default value is {"_x":0, "_y":1, "_z":2}; those are the three Cartesian spatial directions.
    @return Does not return anything. But two binary files are writen instead.
    """
    k_indices, length_dipole_elements = extract_length_dipole_elements(path_to_kfiles, u_to_idx)
    save_pickle(k_indices, os.path.join(output_path, kindx_pkl_fname))
    save_pickle(length_dipole_elements, os.path.join(output_path, lendipel_pkl_fname))

def load_kpoints_data(kdata_fnames={}, output_path=".", kdata_default_fnames={"k_indices":"k_indices.pkl", "length_dipole_elements":"length_dipole_elements.pkl"}):
    """! Loads objects from pickled (binary) files with data (e.g., python data structure or NumPy array) extracted from k-points output files.
    @param kdata_fnames Dictionary (key:value pairs) where values are the paths to the pickled files to be loaded; the loaded object will be stored into the output dictionary with the same corresponding key. Default value is an empty dictionary; if some entries are provided, then those will overwrite matching entries in 'kdata_default_fnames'.
    @param output_path Path to the folder where the output pickled (binary) files should be stored. The default value is the current directory of the function execution. (".").
    @param kdata_default_fnames Dictionary with default names of encrypted files to load as values; regarding keys, it works the same way as 'kdata_fnames', but the latter has priority; default entries will be loaded if available, if not overwritten by matching entries in 'kdata_fnames'. Default value is {"k_indices":"k_indices.pkl", "length_dipole_elements":"length_dipole_elements.pkl"}.
    @return output_kdata Output dictionary storing the loaded objects (using the keys from 'kdata_fnames' and/or 'kdata_default_fnames' (read their descriptions to learn more).
    """
    output_kdata = {}
    kdata_default_fnames.update(kdata_fnames) ; kdata_fnames = kdata_default_fnames
    for kdata_label,kdata_fname in kdata_fnames.items():
        kdata_path = os.path.join(output_path, kdata_fname)
        if not os.path.exists(kdata_path): continue
        output_kdata[kdata_label] = load_pickle(kdata_path)
    return output_kdata

def extract_bandstructure(bandstructure_path):
    """! Extracts useful data (e.g., python data structure or NumPy array) from the output file named 'bandstructure' by default by Octopus; this file contains information related to the 1D k-points path specified in the input file.
    @param Path to the bandstructure output file to be parsed and from which data will be extracted.
    @return Dictionary with the extracted data:
        'k_path' (1: NumPy array detailing the position of each k-point along the 1D path, ranging from 0 -starting point- to 1 -ending point-)
        'n_kpoints' (2: number of k-points along the 1D path)
        'k_coords' (3: reduced coordinates of the k-points along the 1D path)
        'state_energies' (4: energies/eigenvalues of the 1-electron states for each of the k-points along the 1D path)
    """
    extracted_band = {}
    data = np.loadtxt(bandstructure_path)
    # https://www.octopus-code.org/documentation/main/tutorial/periodic_systems/periodic_systems/
    extracted_band["k_path"] = data[:,0]
    extracted_band["n_kpoints"] = data.shape[0]
    extracted_band["k_coords"] = data[:,1:4] # Reduced coordinates
    extracted_band["state_energies"] = data[:,4:]
    return extracted_band

def extract_info(info_path, bandstructure_path=None):
    """! Extracts useful data (e.g., python data structure or NumPy array) from the output file named 'info' by default by Octopus; this file contains the main information and results directly resulting from the Octopus calculation. The extraction can also be combined the parsing of output file named 'bandstructure' by default by Octopus (read the documentation of the 'extract_bandstructure' function to learn more).
    @param info_path Path to the info output file to be parsed and from which data will be extracted.
    @param Path to the output file named 'bandstructure' by default by Octopus; this argument is optional, but if provided, it will change the structure of the output dictionary in order to include a higher level of detail, by differentiating between 3D k-points (key '3D'), the 1D k-points path specified in the input file (key '1D') and the more general information (key 'global'). Default value is None.
    @return Dictionary with the extracted data (regardless of 'bandstructure_path'):
        (Provided in '1D' and '3D' if 'bandstructure_path' is specified:)
        'n_kpoints' (1: number of loaded k-points)
        'k_weights' (2: weight associated with each k-points; apparently set to 0 for the 1D path)
        'k_coords' (3: reduced coordinates of the k-points)
        'occupations' (4: electronic occupation of the 1-electron states, for each k-point; can range from 0 -unoccupied- to 2 -when there is no explicit description of spin)
        'state_energies' (5: energies/eigenvalues of the 1-electron states, for each k-point)
        (Provided in 'global' if 'bandstructure_path' is specified:)
        'n_states' (6: number of 1-electron states)
        'n_ions' (7: number of ions/atoms per unit cell)
        'elements' (8: list of atomic elements for the atoms in the unit cell)
        'ion_forces' (9: forces exerted on each ion in the ground state solution)
        Provided in '1D' only if 'bandstructure_path' is specified:
        'k_path' (10: NumPy array detailing the position of each k-point along the 1D path, ranging from 0 -starting point- to 1 -ending point-)
    """
    extracted_info = {}
    with open(info_path, "r") as infofile: infolines = infofile.readlines()
    for idx,iline in enumerate(infolines):
        # >> Grid information <<
        # ... work in progress ...
        # >> Symmetries <<
        # ... work in progress ...
        # >> Brilloin zone <<
        if "Total number of k-points" in iline:
            extracted_info["n_kpoints"] = int(iline.split("=")[-1])
            continue
        elif "List of k-points:" in iline:
            if "n_kpoints" not in extracted_info: continue
            else: n_kpoints = extracted_info["n_kpoints"]
            from_idx = idx + 3 # (header)
            to_idx = from_idx + n_kpoints
            str_ktable = " ".join(infolines[from_idx:to_idx]).replace("\n","")
            ktable = np.array(str_ktable.split(), dtype=float).reshape((n_kpoints,5))
            extracted_info["k_weights"] = ktable[:,-1]
            extracted_info["k_coords"] = ktable[:,1:-1]
            str_ktable, ktable = None, None # (frees memory)
        # >> Eigenvalues <<
        elif "Eigenvalues [H]" in iline:
            if "n_kpoints" not in extracted_info: continue
            else: n_kpoints = extracted_info["n_kpoints"]
            from_idx = idx + 2 # (header ; #k lines are temp. kept)
            for jdx,jline in enumerate(infolines[from_idx+1:]):
                if jline.strip().startswith("#k"): n_states = jdx ; break
            extracted_info["n_states"] = n_states
            to_idx = from_idx + (1+n_states) * n_kpoints
            eiglines = [line for line in infolines[from_idx:to_idx] if not line.startswith("#k")]
            str_eigtable = " ".join(eiglines).replace("\n","")
            str_eigtable = str_eigtable.replace("--","nan") # (in case of null spin info)
            eigtable = np.array(str_eigtable.split(), dtype=float).reshape((n_kpoints,n_states,4))
            extracted_info["occupations"] = eigtable[:,:,-1]
            extracted_info["state_energies"] = eigtable[:,:,-2]
            #spin = eigtable[:,:,1] # (this will probably need more work if spin info is available)
            str_eigtable, eigtable = None, None # (frees memory)
        # >> Total energy contributions <<
        # ... work in progress ...
        # >> Convergence <<
        # ... work in progress ...
        # >> Forces on ions <<
        elif "Forces on the ions" in iline:
            from_idx = idx + 2 # (header)
            for jdx,jline in enumerate(infolines[from_idx:]):
                if jline.strip().startswith("----"): n_ions = jdx ; break
            extracted_info["n_ions"] = n_ions
            to_idx = from_idx + n_ions
            str_forcetable = " ".join(infolines[from_idx:to_idx]).replace("\n","")
            forcetable = np.array(str_forcetable.split()).reshape((n_ions,5))
            extracted_info["elements"] = forcetable[:,1]
            extracted_info["ion_forces"] = forcetable[:,2:5].astype(float)
            str_forcetable, forcetable = None, None # (frees memory)
    # >> Optional separation between 3D and 1D k-points <<
    if bandstructure_path:
        extracted_band = extract_bandstructure(bandstructure_path)
        n_kpoints_1d = extracted_band["n_kpoints"]
        # 1D k-point path
        extracted_1d = extracted_band
        extracted_1d["k_weights"] = extracted_info["k_weights"][-n_kpoints_1d:]
        extracted_1d["occupations"] = extracted_info["occupations"][-n_kpoints_1d:]
        # 3D k-points
        extracted_3d = {}
        extracted_3d["n_kpoints"] = extracted_info["n_kpoints"] - n_kpoints_1d
        extracted_3d["k_coords"] = extracted_info["k_coords"][:-n_kpoints_1d]
        extracted_3d["state_energies"] = extracted_info["state_energies"][:-n_kpoints_1d]
        extracted_3d["k_weights"] = extracted_info["k_weights"][:-n_kpoints_1d]
        extracted_3d["occupations"] = extracted_info["occupations"][:-n_kpoints_1d]
        # Global info
        extracted_global = {}
        extracted_global["n_states"] = extracted_info["n_states"]
        extracted_global["n_ions"] = extracted_info["n_ions"]
        extracted_global["elements"] = extracted_info["elements"]
        extracted_global["ion_forces"] = extracted_info["ion_forces"]
        # Combined extracted info
        extracted_info = {"1D":extracted_1d, "3D":extracted_3d, "global":extracted_global}
    # >> Returns dictionary with extracted info <<
    return extracted_info

##############################
#   Andrés Ignacio Bertoni   #
# (andresibertoni@gmail.com) #
##############################
