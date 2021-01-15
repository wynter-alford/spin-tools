# spin-tools
Tools for measuring fidelity of a pulse sequence and calculating magnus expansion terms

# Spin Simulation
Simulates a given pulse sequence that is entered in terms of its unitary operators.  The actual sequence unitary after a certain number of cycles (currently 10 milliseconds) will be compared with the expected unitary based on an inputted average Hamiltonian for the sequence to obtain a fidelity score.  There are four input parameters: offset and chemical shift which together form the parameter Delta; the length of the rf pulses; the spacing between the pulses; and the coupling strength between the spins.  The number of spins and grid meshing can also be adjusted, though this increases the runtime considerably (exponentially in the case of the number of spins).

# Magnus
Numerically computes the first five terms (order 0 through 4) of the magnus expansion for a given rf pulse sequence assuming instantaneous pulses, and also computes the magnitude of those terms.
