#!/usr/bin/env python
# -*- coding: utf-8 -*-

## @package octopus_slabs.Libs.libElectrical
# Module with functions to compute and plot the complex electrical conductivity matrix and the complex electrical permittivity, from transition electric dipole moment matrix elements.
# @author Module written by Andrés I. Bertoni (andresibertoni@gmail.com).

import numpy as np

# >>> AUXILIARY functions <<<

def sum_y_same_x(xy_points, round_to_decimals=6):
    """! Converts a set of xy-points, which may be repeated, into a set of points with unique x-values. In the case of repeated x-values, the corresponding y-values are summed.
    @param xy_points Set of xy points, where many instances of the same x-value can be encountered.
    @param round_to_decimals Resolution of x-values (number of decimal places); important for recognizing degeneracy (recognizing values as equal). Default value is 6 decimal places.
    @return A 2D NumPy array with shape (n,2), where n is the number of points with unique x.
    """
    # The points are ordered from the smallest to the largest x-value
    order = np.argsort(xy_points[:,0])
    x_points = np.round(xy_points[order,0], round_to_decimals)
    y_points = xy_points[order,1]
    # y-values are grouped and summed for repeated x-values
    unique_x, firstenc_indices = np.unique(x_points, return_index=True)
    grouped_y = np.split(y_points, firstenc_indices[1:])
    sum_y = np.array([np.sum(y_group) for y_group in grouped_y])
    # xy-points with unique x values are outputed as a NumPy array (with two columns)
    uni_sum_xy_points = np.vstack([unique_x, sum_y]).T
    return uni_sum_xy_points

# >>> CONDUCTIVITY: Kubo-Greenwood Approach <<<

def internal_get_args_conductivity_kg(kpoints_data, info_data):
    """! Extracts the arguments needed for the internal call of the function for computing the conductivity (with the Kubo-Greenwood approach) from the preloaded data structures. Returns them in the correct order so that the calls can be nested in a single line (remember to unpack the tuple for the function call).
    @param kpoints_data Data extracted from k-points output files (e.g., transition electric dipole moment matrix elements). See the 'extract_transition_dipole' function in the 'libOctopusIO' module for extended documentation.
    @param info_data Data extracted from the 'info' output file. See the 'extract_info' function in the 'libOctopusIO' module for extended documentation.
    @return A tuple with the arguments (conveniently ordered) needed for the 'internal_get_conductivity_kg' function.
    """
    n_kpoints = info_data["3D"]["n_kpoints"]
    transition_dipole = kpoints_data["transition_dipole"][:n_kpoints]
    state_energies = info_data["3D"]["state_energies"]
    occupations = info_data["3D"]["occupations"]
    cell_volume = info_data["global"]["cell_volume"]
    k_weights = info_data["3D"]["k_weights"]
    return transition_dipole, state_energies, occupations, cell_volume, k_weights

def internal_get_conductivity_kg(transition_dipole, state_energies, occupations, cell_volume, k_weights=None, damping_parameter=0, round_to_decimals=5, negative_contributions=True):
    """! Computes the complex electrical conductivity matrix with the Kubo-Greenwood approach, from the ground state transition electric dipole moment matrix elements.
    Technical details:
        - Firstly, the conductivity contributions (cc) are computed from the ground state transition electric dipole moment matrix elements, for the full 3D Brillouin k-space, for the frequencies (w) corresponding to transitions from partially occupied states to partially unoccupied states; conductivity contributions for degenerate transitions are summed up, producing a set of (w',cc) points with unique frequencies w'.
        - Secondly, each of the points from the first part, which would correspond each to a Dirac delta for a perfect crystal at 0 K, is replaced by a Lorentzian function with an arbitrary width parameter that we've named 'damping_parameter' as it has physical meaning (adds width to the state energies) for non-perfect crystals and/or non-zero temperature; however, this parameter also serves as a numerical strategy to deal with the fact that only a finite number of Dirac deltas are being computed (because of the basis set truncation, as infinite are required to reproduce the continuous function); the set of points is converted into a numerical function of summed lorentzian conductivity contributions in a grid of frequencies going from zero to the largest transition frequency plus 7 gamma (as +7 gamma includes ~95.5% of the area of the last peak).
        - To deal with the singularity at frequency zero, it is necessary to include the contribution from negative frequencies (this is related to the partity of the real/imaginary parts of the conductivity); doing this also corrects the low-frequency contributions to the conductivity.
        - Formulas (in the velocity/momentum gauge) are available in the following paper: L. Calderín, V.V. Karasiev, S.B. Trickey. "Kubo-Greenwood electrical conductivity formulation and implementation for projector augmented wave datasets". Computer Physics Communications 221, pp. 118-142 (2017). The formulas in this implementation are in the length gauge.
        - WARNING! The formula used for the conductivity is the following: σ(ω) = Re{σ(ω)} - i Im{σ(ω)}. Notice the negative sign; with this convention, electrical conductivity results are consistent with experimental data from Palik tables for bulk Si and Mo; see also http://lampx.tugraz.at/~hadley/memm/materials/silicon/silicon.php for silicon.
    @param transition_dipole Ground state transition electric dipole moment matrix elements, for the full 3D Brillouin k-space. ATTENTION: the transition electric dipole matrix elements should be provided in the length gauge and in atomic units.
    @param state_energies Energies of the 1-electron states for each of the k-points of the full 3D Brillouin k-space.
    @param occupations Electronic occupation of the 1-electron states, for each of the k-points of the full 3D Brillouin k-space; can range from 0 -unoccupied- to 2 -when there is no explicit description of spin-.
    @param cell_volume Volume of the unit cell in atomic units.  
    @param k_weights Weight associated with each of the k-points of the full 3D Brillouin k-space; used for integrating the full 3D Brillouin k-space. Default value is None; in that case it is assumed that all k-points have the same weight (with the sum of weights equal to unity).
    @param damping_parameter Width of the Lorentzian functions that approximate the Dirac deltas, in atomic units of energy; needed both for numerical representation of the conductivity (under the basis set trunctation approximation) and in the case of non-perfect crystals and/or non-zero temperature. It is recommended to compute the conductivity using several values for this parameter (e.g. between the equivalent to 0.1 and 1 eV). Default value is 0.
    @param round_to_decimals Resolution of the transition frequencies (important for recognizing degeneration) and also for the grid; if the resolution is two small, the size of the grid may become impractically big and the computation time will be unnecessarily longer. Default value is 5 decimal places (resolution of 0.00001 or 1e-5).
    @param negative_contributions If the contributions from negative frequencies should be included to correct the conductivity for lower frequencies and avoid the singularity at zero frequency. Default value is True. We recommend including it where possible, but the switch was included for completeness.
    @return Returns both a grid of frequencies of an hypothetical external field (in atomic units of energy) and the complex electrical conductivity matrix as computed by the Kubo-Greenwood approach; it is a 3x3 (αxβ) matrix that encodes the current response of the system along a certain spatial direction (α) due to an external field perturbation along an arbitrary spatial direction (β); each element of the conductivity matrix has the same size as the frequency grid, being a complex function of the frequency. The non-diagonal elements may be erroneous if the ground state transition electric dipole moment matrix elements do not include sign information (artificial positive values).
    """
    # <<< FIRST PART : Points >>>
    # Extracts n_states and n_kpoints if they were not provided
    n_kpoints, n_dim, n_states = transition_dipole.shape[:-1]
    # Checks for the k-point integration weights array
    # (if all weights are the same, then it is better to compute this value by hand -more precise-)
    nonexistent_weights = (k_weights is None)
    equal_weights = (np.all(k_weights == k_weights[0]))
    wrong_weights = (np.sum(k_weights) != 1)
    if nonexistent_weights or equal_weights or wrong_weights:
        k_weights = np.ones(n_kpoints) * (1/n_kpoints)
    # Constant prefactor
    prefactor = 2 / cell_volume
    # The prefactor should probably include a 1/4 ((1/2)²) factor to fix issue with Octopus implementation
    #(1/2) factor missing during the conversion from the gradient/velocity gauge to the length gauge,
    #  for each transition dipole element...)
    conductivity = []
    # (At the end, the conductivity will be a (n_dim x n_dim) matrix made of nested lists,
    # with NumPy arrays as elements; each array is a set of (omega_mn, conductivity) points.) 
    for u_idx in range(n_dim):
        conductivity.append([])
        for u_jdx in range(n_dim):
            conductivity[u_idx].append([])
            for kdx,w_k in enumerate(k_weights):
                # Occupation difference : Δf_nm = f(En) - f(Em)
                delta_occ_nm = np.subtract.outer(occupations[kdx], occupations[kdx])
                allowed = tuple(map(tuple, np.argwhere(delta_occ_nm > 0).T))
                delta_occ_nm = delta_occ_nm[allowed]
                # Energy difference : w_nm = (Em - En)/ℏ
                omega_mn = np.subtract.outer(state_energies[kdx], state_energies[kdx]).T
                omega_mn = omega_mn[allowed]
                # k-point integration weight
                w_k = k_weights[kdx]
                # Transition dipole momenta
                trdip_idx_mn = transition_dipole[kdx,u_idx].T # <m|µi|n>
                trdip_jdx_nm = transition_dipole[kdx,u_jdx]   # <n|µj|m>
                trdip_idx_mn = trdip_idx_mn[allowed]
                trdip_jdx_nm = trdip_jdx_nm[allowed]
                # Computing conductivity points (because summands are Dirac deltas)
                conductivity_values = prefactor * w_k * delta_occ_nm * trdip_idx_mn * trdip_jdx_nm * omega_mn
                conductivity_points = zip(omega_mn, conductivity_values)
                conductivity[u_idx][u_jdx].extend(conductivity_points)
            # Converting conductivity_ij point-sets into numpy arrays
            conductivity[u_idx][u_jdx] = np.array(conductivity[u_idx][u_jdx])
            # Sum conductivities with the same transition frequency
            conductivity[u_idx][u_jdx] = sum_y_same_x(conductivity[u_idx][u_jdx], round_to_decimals)
    # <<< SECOND PART : Grids >>>
    # The information contained in the set of points is transfered to a numerical grid
    # Starts by creating a common grid for omega
    max_omega = np.max([np.max(conductivity[u_idx][u_jdx][:,0]) for u_idx in range(n_dim) for u_jdx in range(n_dim)]) + (damping_parameter * 7) # (±7 gamma includes ~91% of the area of a Lorentzian peak) 
    resolution = 1/(10**round_to_decimals)
    grid_points = 1 + int(max_omega/resolution)
    omega_grid = np.round(np.linspace(0, max_omega, grid_points), round_to_decimals)
    # Creates and loads the complex arrays for the elements of the conductivity
    # (overwrites pre-existing nested-lists to save memory)
    for u_idx in range(n_dim):
        for u_jdx in range(n_dim):
            omega_ij, conduct_ij = conductivity[u_idx][u_jdx].T
            conduct_grid_ij = np.zeros(omega_grid.size, dtype=complex)
            # Without damping, real part is a set of Dirac deltas
            if damping_parameter == 0:
                point_indices = np.argwhere(np.in1d(omega_grid, omega_ij)).flatten()
                conduct_grid_ij[point_indices] += conduct_ij / resolution # WARNING! For some reason that we still have not been able to recognize, this requires dividing by the energy resolution so these results coincide with the Lorentzian approach...
                iterate_real = False
            else: iterate_real = True
            # Using (broken: δ→0) lorentzian-like functions with(out) damping
            # (this part is slow but inexpensive memory-wise)
            for kdx,omega_ijk in enumerate(omega_ij):
                real_num = (damping_parameter/2)
                imag_num_i = (omega_ijk - omega_grid)
                denominator_i = (real_num**2) + (imag_num_i**2)
                if negative_contributions:
                    neg_imag_num_i = (omega_ijk + omega_grid)
                    neg_denominator_i = (real_num**2) + (neg_imag_num_i**2)
                # Real part
                if iterate_real:
                    real_part_i = conduct_ij[kdx] * real_num/denominator_i
                    if negative_contributions:
                        real_part_i += - conduct_ij[kdx] * real_num/neg_denominator_i
                        real_part_i *= 1/2
                    conduct_grid_ij += real_part_i
                # Imag part
                imag_part_i = conduct_ij[kdx] * imag_num_i/denominator_i
                if negative_contributions:
                    imag_part_i += - conduct_ij[kdx] * neg_imag_num_i/neg_denominator_i
                    imag_part_i *= 1/2
                if damping_parameter == 0:
                    imag_part_i = np.where(np.isnan(imag_part_i), 0, imag_part_i)
                conduct_grid_ij += - 1j * imag_part_i # WARNING: Negative sign present!
                # Check for Silicon: http://lampx.tugraz.at/~hadley/memm/materials/silicon/silicon.php
            # Stores the ij element of the complex conductivity
            conductivity[u_idx][u_jdx] = conduct_grid_ij
    # Returns the complex conductivity as a single NumPy matrix (alongside the omega grid)
    conductivity = np.array(conductivity, dtype=complex)
    return omega_grid, conductivity

# >>> CONDUCTIVITY: Drude Approach <<<

def internal_get_conductivity_d(omega_grid, eff_charge_to_mass, damping_parameter, round_to_decimals=5):
    """! Computes the complex electrical conductivity matrix with the Drude approximation. This approach requires an inputed parameter (which is usually obtained by fitting to experimental data).
    Technical details:
        The formula for the Drude expression (derived from the complex conductivity can be found in the following paper: L. Calderín, V.V. Karasiev, S.B. Trickey. "Kubo-Greenwood electrical conductivity formulation and implementation for projector augmented wave datasets". Computer Physics Communications 221, pp. 118-142 (2017).
        Calderín et al. only show the Drude contribution for the real part of the complex conductivity; the imaginary part may come from the Onsager reciprocal relations.
    @param omega_grid Corresponds to the grid of frequencies of an hypothetical external field (in atomic units of energy). This can be obtained as the output from the calculation of the conductivity matrix using another approach, such as the Kubo-Greenwood approach; this is important if conductivity contributions are to be easily summed.
    @param eff_charge_to_mass The effective charge/mass ratio is the main parameter of the Drude expression, which is expected to capture contributions to the complex conductivity that are given by (i) the degenerate 'vertical' inter-band transition and (ii) the intra-band transitions. As far as we know, these contributions cannot be extracted from non-TD calculations for obtaining the ground state transition electric dipole moment matrix elements. This implementation allows this quantity to be a complex number, leaving aside its physical significance.
    @param damping_parameter Width of the Lorentzian functions that approximate the Dirac deltas, in atomic units of energy; needed both for numerical representation of the conductivity (under the basis set trunctation approximation) and in the case of non-perfect crystals and/or non-zero temperature. It is recommended to compute the conductivity using several values for this parameter (e.g. between the equivalent to 0.1 and 1 eV).
    @param round_to_decimals Resolution of the transition frequencies (important for recognizing degeneration) and also for the grid; if the resolution is two small, the size of the grid may become impractically big and the computation time will be unnecessarily longer. Default value is 5 decimal places (resolution of 0.00001 or 1e-5).
    @return The complex electrical conductivity matrix as calculated by the Drude approximation. The diagonal elements are all equal and the non-diagonal elements are set to zero. It is returned as a 3x3 (αxβ) matrix, since it is meant to encode the current response of the system along a given spatial direction (α) to an external field perturbation along an arbitrary spatial direction (β); each element of the conductivity matrix has the same size as the frequency grid, being a complex function of frequency.
    """
    # Checks the validity of the damping parameter
    if damping_parameter == 0: raise RuntimeError("The damping parameter cannot be zero! It should be twice the average inter-collision frequency in the context of the Drude model (i.e., the reciprocal of the average inter-collision time).")
    # Computes the Drude conductivity as a complex 3x3 matrix
    tau = 2/damping_parameter # (tau is the average inter-collision time)
    n_points = len(omega_grid)
    conductivity = np.zeros((3,3,n_points), dtype=complex)
    for idx in range(3):
        conductivity[idx,idx,:] = np.ones(n_points)
    if type(eff_charge_to_mass) != complex:
        # (it is converted into a complex number with equal real and imaginary parts)
        eff_charge_to_mass *= (1+1j)
    numerator = eff_charge_to_mass * tau
    denominator = 1 + ((tau * omega_grid)**2)
    conductivity *= (numerator/denominator)
    conductivity = np.round(conductivity, round_to_decimals)
    # Returns the complex conductivity
    return conductivity

# >>> MAIN functions <<<

def get_conductivity(data, mode="Kubo-Greenwood", damping_parameter=0, round_to_decimals=5, negative_contributions=True):
    """! Computes the complex electrical conductivity matrix from the ground state transition electric dipole moment matrix elements. Currently the only implemented mode is the Kubo-Greenwood approach.
    @param data This is a dictionary with the data needed for the computation.
    Required entries vary with the selected mode:
    -> Kubo-Greenwood:
        "kpoints" : Data extracted from k-points output files (e.g., transition electric dipole moment matrix elements). See the 'extract_transition_dipole' function in the 'libOctopusIO' module for extended documentation.
        "info" : Data extracted from the 'info' output file. See the 'extract_info' function in the 'libOctopusIO' module for extended documentation.
    -> Drude:
        "omega_grid" : Grid of frequencies that can be obtained as the output of the computation of the conductivity matrix with another approach, such as the Kubo-Greenwood approach; this is important if conductivity contributions are to be easily summed.
        "eff_charge_to_mass" : The effective charge/mass ratio is the main parameter of the Drude expression (see 'internal_get_conductivity_d' for extended documentation).
    @param mode How to compute the complex electrical conductivity. Default value is "Kubo-Greenwood" as it is currently the only mode implemented (other valid values are "kg", capital letters, spaces and hyphens are ignored). Read the 'internal_get_conductivity_kg' function for extended documentation.
    @param damping_parameter Parameter needed both for numerical representation of the conductivity and in the case of non-perfect crystals and/or non-zero temperature. Read the 'internal_get_conductivity_kg' function for extended documentation. Default value is 0.
    @param round_to_decimals Frequency resolution. Default value is 5 decimal places (resolution of 0.00001 or 1e-5). Note that very small resolutions considerably increase the computation time and memory required, without necessarily 'visually' improving the quality of the results. Read the function 'internal_get_conductivity_kg' for extended documentation.
    @param negative_contributions Whether the conductivity should be corrected for lower frequencies and to avoid the singularity at zero frequency. Default value is True; do not change if you are unsure. Read the function 'internal_get_conductivity_kg' for extended documentation. (This is not needed for the Drude conductivity!)
    @return Returns both the frequency grid (in atomic units of energy) and the complex electrical conductivity matrix (3x3, also in atomic units), that encodes the current response of the system due to an external field perturbation; each element of the conductivity matrix has the same size as the frequency grid, being a complex function of the frequency. Read the function 'internal_get_conductivity_kg' for extended documentation.
    """
    cmode = mode.strip().lower().replace(" ","").replace("-","")
    if cmode in ["kubogreenwood", "kg"]:
        kpoints_data = data["kpoints"]
        info_data = data["info"]
        omega_grid, conductivity = internal_get_conductivity_kg(*internal_get_args_conductivity_kg(kpoints_data, info_data), damping_parameter, round_to_decimals, negative_contributions)
    elif cmode in ["drude", "d"]:
        omega_grid = data["omega_grid"]
        eff_charge_to_mass = data["eff_charge_to_mass"]
        conductivity = internal_get_conductivity_d(omega_grid, eff_charge_to_mass, damping_parameter, round_to_decimals)
    return omega_grid, conductivity

def get_permittivity(omega, conductivity, diagonal_average=True):
    """! Computes the complex electrical permittivity from the complex electrical conductivity.
    Employs the folowing formula:
        ε(ω) = 1 - σ(ω)/(i.ω.ε0) --> ε(ω) = 1 + i . σ(ω)/(ω.ε0) , where ε0 = 1/(4π)
    Using this sign convention and the sign convention for the real and imaginary parts of the conductivity (σ(ω) = Re{σ(ω)} - i Im{σ(ω)}, see 'internal_get_conductivity_kg' for extended documentation), electrical permittivity results are consistent with experimental data from Palik tables for bulk Si and Mo.
    @param omega Grid of frequencies in atomic units of energy.
    @param conductivity Complex electrical conductivity matrix (3x3) in atomic units.
    @param diagonal_average Whether isotrophy is assumed and, therefore, the average of the diagonal elements of the complex electrical conductivity is used. Default value is True (isotrophy is assumed).
    @return Returns the complex electrical permittivity.
    """
    epsilon = 1 + 1j * (conductivity/omega) * (4*np.pi)
    if diagonal_average: epsilon = np.trace(epsilon)/3
    return epsilon

# <<< EXTRA Drude function (for fitting!) >>>

from scipy.optimize import curve_fit

def fit_drude_to_permittivity(omega_grid, base_conductivity, target_permittivity, damping_parameter, fit_omega_range=None, target_param_bounds=(0,1), n_inits=20, round_to_decimals=5, verbose=True):
    """! Finds the Drude parameter (i.e., effective charge/mass) that produces the Drude contribution to the complex electrical conductivity that, together with the inputed base conductivity, best fits the also given target electrical permittivity.
       Technical details:
           - The main formula is available in the following paper: L. Calderín, V.V. Karasiev, S.B. Trickey. "Kubo-Greenwood electrical conductivity formulation and implementation for projector augmented wave datasets". Computer Physics Communications 221, pp. 118-142 (2017).
           - The following diagram may be of help:
                 (Drude parameter --> Drude conductivity) + Base conductivity (...)
                 (...) --> Base+Drude permittivity <--(comparison)--> Target permittivity
    @param omega_grid Corresponds to the grid of frequencies of an hypothetical external field (in atomic units of energy). This array should be the same size as the specified target permittivity and the last dimension of the base conductivity.
    @param base_conductivity This is the base conductivity to be added to the calculated Drude contribution to the complex electrical conductivity before conversion to an electrical permittivity that will be compared to the given target; can be obtained as the output from the calculation of the conductivity matrix using an approach such as the Kubo-Greenwood one.
    @param target_permittivity This is the target permittivity to which the Drude contribution is fitted. This is usually data interpolated from experimental tables.
    @param damping_parameter Width of the Lorentzian functions approximating the Dirac deltas, in atomic units of energy; needed for the numerical representation of the conductivity and for the estimation of the average intercollision time in the Drude model.
	@param fit_omega_range If specified, it will use only the specified range (in omega units) to perform the fit. It defaults to None, which means that the entire domain will be used in the fitting procedure.
    @param target_param_bounds Bounds for the Drude parameter to be optimised (unitless). It defaults to (0,1), which means that the minimum possible value of the optimised parameter is 0 and its maximum is 1; by default, this parameter is always positive.
    @param n_inits Number (integer) of optimisations to be performed from different initial guesses. These guesses are extracted from the specified bounds of the target parameter (excluding the extremal values). This number should be considerably large to ensure a successful optimisation; however, the total computational time will increase linearly. The default value is 20, but a value of 100 would be ideal (although more expensive).
    @param round_to_decimals Resolution of the transition frequencies and the output fitted parameter; if the resolution is two small, the size of the grid may become impractically big and the computation time will be unnecessarily longer. Default value is 5 decimal places (resolution of 0.00001 or 1e-5).
    @param verbose Whether or not to print messages at runtime. The default is True, i.e. all messages are displayed by default.
    @return Returns the optimum Drude parameter (i.e. the effective charge to mass ratio), which can be used, for example, to calculate the Drude contribution to the complex electrical conductivity.
    """
    # Truncates fitting region
    if fit_omega_range:
        from_idx = np.argmin(np.abs(omega_grid - fit_omega_range[0]))
        to_idx = np.argmin(np.abs(omega_grid - fit_omega_range[1]))
    else: from_idx, to_idx = None,None
    tr_x = omega_grid[from_idx:to_idx]
    tr_y_target = target_permittivity[from_idx:to_idx]
    # Defines the model permittivity (from conductivity)
    def model_permittivity(view_omega_grid, eff_charge_to_mass):
        omega_grid = view_omega_grid[:view_omega_grid.size//2]
        # Following arguments do not enter explicitly:
        #   base_conductivity, from_idx, to_idx
        #   damping_parameter, round_to_decimals
        # (sigma : conductivity | epsilon : permittivity
        base_sigma = base_conductivity[:,:,from_idx:to_idx]
        data = {"omega_grid":omega_grid, "eff_charge_to_mass":eff_charge_to_mass}
        drude_sigma = get_conductivity(data, "Drude", damping_parameter, round_to_decimals)[-1]
        model_sigma = base_sigma + drude_sigma
        model_epsilon = get_permittivity(omega_grid, model_sigma)
        return model_epsilon.view(float)
    # Finds the optimal Drude parameter (i.e. effective charge/mass)
    view_tr_x = (tr_x + 1j * tr_x).view(float)
    view_tr_y_target = tr_y_target.view(float)
    # (tests for multiple initial parameters)
    all_opt_params = [] ; all_pcov_params = []
    for init_param in np.linspace(*target_param_bounds, n_inits+2)[1:-1]:
        fit_params, pcov = curve_fit(model_permittivity, view_tr_x, view_tr_y_target, p0=init_param, bounds=target_param_bounds, method="trf")
        if pcov == 0: continue
        all_opt_params.append(fit_params[0])
        all_pcov_params.append(pcov.flatten()[0])
    opt_param = all_opt_params[np.argmin(all_pcov_params)]
    # Returns the optimal parameter
    if verbose: print(f"Optimal eff. charge/mass: {opt_param:.{round_to_decimals}f}")
    return np.round(opt_param, round_to_decimals)

# <<< EXTRA general functions >>>

import matplotlib.pyplot as plt
import scipy.constants as spc

def convert_eV_to_nm(energy_eV):
    """! Converts energy in electron-volt units (eV) to wavelength in units of nanometers (nm). It works with arrays and in the opposite direction of conversion (from nm to eV) without having to specify it.
    @param energy_eV Single value or array of energy values in eV units to be converted to wavelength in nm. If wavelengths in nm are provided instead, these are converted to energies in eV.
    @return Returns the converted value or array (read above).
    """
    c = spc.physical_constants["speed of light in vacuum"][0] # m/s
    h = spc.physical_constants["Planck constant"][0] # J.s
    J_to_eV = spc.physical_constants["joule-electron volt relationship"][0] # eV/J
    m_to_nm = 1e9 # nm/m
    hc_eVnm = (c * h * m_to_nm * J_to_eV) # eV.nm
    wavelenght_nm = hc_eVnm / energy_eV # nm
    return wavelenght_nm

def load_epsilon_table(path_to_file, in_angstrom=True, out_energy=True):
    """! From experimental tables containing the wavelength (λ) dependence of the refractive index (n) and the extinction coefficient (k), computes the complex electrical permittivity (ε) as a function of wavelenth in nanometers (nm) or energy in electron-volt (eV) units. A table with three columns (λ | n | k) is expected; the table wavelength can be accepted in nanometers (nm) or in Ångström (Å).
    @param path_to_file Path to the text file with the table (lines commented with '#' are ignored).
    @param in_angstrom Whether the wavelength is in Å; otherwise assumed to be in nm. Default value is True (wavelengths are expected in Å).
    @param out_energy Whether the output grid should be wavelength in nm or energy in eV. Default value is True (output grid is energy in eV).
    @return Returns the grid of wavelength (energy) in nm (in eV) and the complex electrical permittivity extracted from the experimental data table.
    """
    table = np.loadtxt(path_to_file)
    wavelength = table[:,0]
    if in_angstrom:
        A_to_nm = 1/10
        wavelength *= A_to_nm
    n = table[:,1] ; k = table[:,2]
    epsilon = (n + 1j * k)**2
    if out_energy:
        energy = convert_eV_to_nm(wavelength)
        return energy, epsilon
    else: return wavelength, epsilon

def interpolate_permittivity(target_omega, omega, permittivity):
    """! Interpolates an electrical permittivity curve for a denser frequency domain. This is useful for converting data from experimental tables to arrays of sizes that match those coming from calculations (allowing operations such as addition and subtraction, etc.).
    @param target_omega Corresponds to the denser grid of frequencies of an hypothetical external field (in atomic units of energy). This can be obtained as the output from the calculation of the conductivity matrix using an approach such as the Kubo-Greenwood one.
    @param omega This corresponds to the less dense range of frequencies; it is the domain (x) of the curve to be interpolated.
    @param permittivity This corresponds to the less dense array of the complex permittivity; it is the image (y) of the curve to be interpolated.
    @return Returns the interpolated (denser array) complex electrical permittivity.
    """
    x = target_omega
    xp = omega[::-1]
    re_yp = np.real(permittivity)[::-1]
    im_yp = np.imag(permittivity)[::-1]
    y = np.interp(x, xp, re_yp) + 1j * np.interp(x, xp, im_yp)
    return y

# <<< EXTRA Plot functions >>>

def plot_complex_array(grid, complex_array, max_grid=None, stem_alpha=0.2, x_label="x", y_label="y", plt_obj=None):
    """! Generic function to easily plot the real and imaginary parts of a complex array.
    @param grid Array with the domain of the complex function; it has to have the same size as the complex array.
    @param complex_array Complex array to be plotted; it has to have the same size as the grid.
    @param max_grid The function is plotted up to this value (same units as the grid); if this value is None or greater than the maximum value of the grid, the entire function is then be plotted. Default value is None (the entire domain is plotted).
    @param stem_alpha Transparency of the stems in the plot (from 0 -invisible- to 1 -solid-). Default value is 0.2 (20% solid or 80% transparent).
    @param x_label Label for the grid axis (horizontal axis). Default value is "x".
    @param y_label Label for the axis of the function values (vertical axis). Default value is "y".
    @param plt_obj Plot object for setting parameters; useful when working with subplots. Default value is None, which means that plt.gca() is used to set the plot parameters.
    @return It returns nothing; a plot is displayed instead.
    """
    if plt_obj is None: plt_obj = plt.gca()
    x = grid ; y = complex_array
    y_re = np.real(y) ; y_im = np.imag(y)
    colors = ["black", "red"]
    if max_grid:
        xlim_mask = (x <= max_grid)
        x = x[xlim_mask] ; y_re = y_re[xlim_mask] ; y_im = y_im[xlim_mask]
    marker, stem, base = plt_obj.stem(x, y_re, markerfmt=" ", basefmt=" ", label=f"$\mathfrak{{Re}} ({y_label})$") ; plt.setp(stem, "color", colors[0]) ; plt.setp(stem, "alpha", stem_alpha)
    marker, stem, base = plt_obj.stem(x, y_im, markerfmt=" ", basefmt=" ", label=f"$\mathfrak{{Im}} ({y_label})$") ; plt.setp(stem, "color", colors[1]) ; plt.setp(stem, "alpha", stem_alpha)
    plt_obj.set_xlabel(x_label)
    leg = plt_obj.legend(fontsize=16, frameon=False, handlelength=0, labelcolor=colors)
    [leg.legend_handles[p].set_visible(False) for p in range(2)]

def plot_square_complex_matrix(grid, square_complex_matrix, max_grid=None, udx_to_u=["x","y","z"], diagonal_only=False, subplot_figsize=3, store_path=None, stem_alpha=0.2, x_label="x", y_label="y"):
    """! Generic function to easily plot the real and imaginary parts of a complex square matrix; treats each element as a complex array.
    @param grid Array with the domain of the elements of the complex matrix; it has to have the same size as the elements.
    @param square_complex_matrix Complex square matrix to be plotted; each element has to have the same size as the grid
    @param max_grid The function is plotted up to this value (same units as the grid); if this value is None or greater than the maximum value of the grid, the entire function is then be plotted. Default value is None (the entire domain is plotted).
    @param udx_to_u Labels for the N dimensions of the NxN matrix (in order of occurrence) to be included as subscripts in the subplots' legend entries (real and imaginary parts of the complex function). Default value is ["x", "y", "z"]; therefore, the diagonal elements would be "xx", "yy" and "zz"; while the non-diagonal elements would be "xy", "yz", "zx" (and the labels in reverse order).
    @param diagonal_only Whether to display only diagonal elements or all elements of the complex matrix. Default value is False (all elements are shown).
    @param subplot_figsize Each subplot has a size of M by M, where M is this number (in default MatPlotLib units). Default value is 3 (each subplot is 3x3).
    @param store_path Path to store the plot as a file; the desired file extension must be included (e.g. '.jpg', '.png', '.svg', etc.). Default value is None, which means that the plot is only displayed and not stored.
    @param stem_alpha Transparency of the stems in the plot (from 0 -invisible- to 1 -solid-). Default value is 0.2 (20% solid or 80% transparent).
    @param x_label Label for the grid axis (horizontal axis). Default value is "x".
    @param y_label Label for the subplots' legend entries (real and imaginary parts of the complex function). Default value is "y".
    @return It returns nothing; a plot is displayed instead.
    """
    n_dim = len(udx_to_u)
    if diagonal_only: fig, axs = plt.subplots(1, n_dim, figsize=(n_dim * subplot_figsize, subplot_figsize), sharey=True)
    else: fig, axs = plt.subplots(n_dim, n_dim, figsize=(n_dim * subplot_figsize, n_dim * subplot_figsize), sharey=True)
    for i in range(n_dim):
        for j in range(n_dim):
            if diagonal_only:
                if i != j: continue
                else: ax_ij = axs[i]
            else: ax_ij = axs[i,j]
            y_label_ij = f"{y_label}_{{{udx_to_u[i]}{udx_to_u[j]}}}"
            plot_complex_array(grid, square_complex_matrix[i][j], max_grid, stem_alpha, x_label, y_label_ij, ax_ij)
    fig.tight_layout()
    if store_path is not None: plt.savefig(store_path)

def plot_conductivity(omega, conductivity, max_omega_eV=20, udx_to_u=["x", "y", "z"], diagonal_only=False, subplot_figsize=3, store_path=None, stem_alpha=0.2, x_label=None, y_label=None):
    """! Plots the real and imaginary parts of the complex electrical conductivity matrix (in atomic units) versus the frequency of an hypothetical external field in electron-volt (eV) units of energy.
    @param omega Grid of frequencies in atomic units of energy. Automatically converted to electron-volts (ev) units.
    @param conductivity Complex electrical conductivity matrix (3x3) in atomic units.
    @param max_omega_eV The function is plotted up to this value of energy in eV; if this value is None or greater than the maximum value of the grid, the entire function is then be plotted. Default value is 20 (each element of the complex electrical conductivity is plotted from 0 eV to 20 eV).
    @param udx_to_u Labels for the N dimensions of the NxN matrix (in order of occurrence) to be included as subscripts in the subplots' legend entries (real and imaginary parts of the complex electrical conductivity matrix). Default value is ["x", "y", "z"]; therefore, the diagonal elements would be "xx", "yy" and "zz"; while the non-diagonal elements would be "xy", "yz", "zx" (and the labels in reverse order).
    @param diagonal_only Whether to display only diagonal elements or all elements of the complex electrical conductivity matrix. Default value is False (all elements are shown).
    @param subplot_figsize Each subplot has a size of M by M, where M is this number (in default MatPlotLib units). Default value is 3 (each subplot is 3x3).
    @param store_path Path to store the plot as a file; the desired file extension must be included (e.g. '.jpg', '.png', '.svg', etc.). Default value is None, which means that the plot is only displayed and not stored.
    @param stem_alpha Transparency of the stems in the plot (from 0 -invisible- to 1 -solid-). Default value is 0.2 (20% solid or 80% transparent).
    @param x_label Label for the grid axis (horizontal axis). Default value is None, meaning that 'ω [eV]' will be used for the horizontal axis label.
    @param y_label Label for the subplots' legend entries (real and imaginary parts of the complex electrical conductivity matrix). Default value is None, meaning that 'σ' will be used as label.
    @return It returns nothing; a plot is displayed instead.
    """
    Eh_to_eV = spc.physical_constants["hartree-electron volt relationship"][0]
    if x_label is None: x_label = "ω [eV]"
    if y_label is None: y_label = "σ"
    plot_square_complex_matrix(omega*Eh_to_eV, conductivity, max_omega_eV, udx_to_u, diagonal_only, subplot_figsize, store_path, stem_alpha, x_label, y_label)

##############################
#   Andrés Ignacio Bertoni   #
# (andresibertoni@gmail.com) #
##############################
