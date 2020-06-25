# Present status of running works

## Ground state preparation

June 22nd 2020

* Band structure was checked. 
* Convergence was well obtained. 
* Fermi energy was obtained. Warning: Tancogne-Dejean indicated that a bug may have been found. A restart can help to check that. No consequences other than for display. 

## Linear optical response of Au bulk

June 22nd 2020

* Obtained very nice optical response of Au sample if using the oct-dielectric tool from Octopus. Results are well comparable to experimental results from Johnson and Christie. However, the optical response in IR and mid-IR is not satisfactory. It is unclear if DFT could provide better. 
* Home made routine based on Ohm law has been developed but fails epically. Sato employs Ohm law with a mask function, Yabana employs another method with vector potential A. TO BE CLARIFIED TO SATISFY PROJECT IT4I and ATLANTIC WP5.

## Non-linear optical response for Au bulk
June 22nd 2020
* Procedure for obtaining optical response is dependent to the one of linear response.

## Energy absorption by Au in high intensity regime

June 22nd 2020

* Principle: Experimental observation of reduction of reflectivity must correspond to an increase of absorptivity at high fluence.

### 1. Express laser energy we put in there. 

June 22nd 2020.

* Problem is to express energy we put into the Hamiltonian. Which term is affected into the Hamiltonian? Octopus code is written in a very implicit manner (difficult to find where it could be). Boundary condition conditions of Maxwell equations (given by continuity of macroscopic fields) provide formula given by Yabana PRB 2012. Generalizing this formula is difficult due to temporal integration of J(t). 

### 2. Temporal integration of J(t) for some runs. 

June 22nd 2020. 

* Restart some new runs with outputing J(t) more often. Convergence of P(t) should be obtained (with two methods). 
* If fails, then prepare runs without induced fields => field definition will be more clear.
* If fails, one would have to identify which term we introduce into Hamiltonian, directly looking at Octopus.  

### 3. The role of dressing

June 22nd 2020.

* Dressing must be very important (according to Petek Nat. Comm. 2020). Compute quasi-band structure of Au. Understand the origin of threshold for Energy_eigenvalues getting suddenly from negative to positive values. 
