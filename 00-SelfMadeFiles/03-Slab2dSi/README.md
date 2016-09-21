# Sample

In this folder, sample is a 2D-periodic Si slab sample. 

* Periodic boundary conditions are applied in X and Y directions. 

* Absorbing boundary condition are applied in the Z direction.

# Expected results

Outcome of simulations with this sample are the following. 

* Attempt to reproduce optical data from Palik for Si <100> using two methods: 

1. TD-DFT: this one gives error message. "Use GaugeVectorField".

1. Casidi equation.

# Problem

Similarly to the 3D-periodic case, we obtain: 
```bash
**************************** FATAL ERROR *****************************
*** Fatal Error (description follows)
*--------------------------------------------------------------------
* From node =    0
*--------------------------------------------------------------------
* Kick cannot be applied in a periodic direction. Use GaugeVectorField instead.
**********************************************************************

```

Hence, we shall rather build a big supercell with many atoms and absobring boundary. 

# Later ? 

We could consider to develop GaugeVectorField with Dr Heiko Appel if necessary. 
