# spin-tools
Tools for measuring fidelity of a pulse sequence and calculating magnus expansion terms.  The most important pieces of code are *fullMagnusOneParam.m* and *oldFullAnalysis.m*.

# Magnus Tools.ipynb
The first section of this notebook outlines the code used to compute the terms of the magnus expansion for various pulse sequences.  The second section contains all of the results generated by simulations so far.

# fullMagnusOneParam.m
This is the actual script used to compute terms in the magnus expansion.

# spinSimOneParam.ipynb
This script generates 1D plots of fidelity from a single parameter for various pulse sequences.

# oldFullAnalysis.m
This is the script used during the summer/fall of 2020 to compute fidelity varying all four parameters.

# fullAnalysis.m
This was intended to work as an updated version of *oldFullAnalysis.m* that was more efficient.  The code to generate the first and second order terms of the average Hamiltonian was also included so that fidelity could be computed using the first three terms to approximate the Hamiltonian instead of just the zeroth order term.  However, the results this code produced when considering fidelity w.r.t. the zeroth order average Hamiltonian did not match up to the results from *oldFullAnalysis.m*.  This discrepancy has not yet been resolved.

# fullMagnusGrid.m
Decrepit.
