# On-the-fly Calculation of Transport Properties of Fluids using the order-n algorithm in Equilibrium Molecular Dynamics


These points should be considered for the correct implementation of fix ordern combined with compute commands:

1) fix ordern always accepts vector as an input from compute commands. This vector consists of 6 components of the stress tensor for viscosity, 6 compnents of the heat flux vector for thermal conductivity, and 5N components of positions, global-ID, and ... of all atoms for diffusion coefficients.

2) only 1 transport property can be computed per fix ordern. Therefore, for computing the diffusion coefficients, viscosities, and thermal conductivity of a system, 3 compute commands and 3 fix order commands are required.

Remained issues: (23-10-2018)

1) Add warning/error handling for all possible situations:
    a) insufficient number of input arguments
    b) incorrect input arguments (e.g. local vs. global)
    c) insufficient memory to store data

2) Modifying the output functions for all transport properties
    a) possible change in the name of the output file if not specified
    b) make a separate function

3) Make use of dynamic arrays for large arrays
    a) possible change of the size during a simulation
    b) specifying the number of buffer sizes for each block
    c) no predefined value for the maximum length of each array

4) Remove all irrelevant attributes/methods/parts
    a) improving the performance by changing unnecessary operations

5) Solve the problem of restart_timestep