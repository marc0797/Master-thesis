# Experiment Reproduction Instructions

To reproduce the experiments in this folder, please follow the steps below:

## Matlab Notebooks

1. Open the `Matlab Notebooks` folder.
2. Ensure that you have MATLAB installed on your system.
3. Run the MATLAB notebooks in the specified order, as mentioned in the experiment documentation.
4. Make sure to provide the required input data or modify the notebook parameters as needed.
5. Save the results and any generated output files in the appropriate location.

### Transverse Field Ising model 1D, 2 spins.

For the first experiment, the Matlab file `Ising1D_Neq2.m` showcases how to
optimize the variational parameters by imposing that the gradients of the energy with respect to
them must be zero.

The solution is found analytically using symbolic variables and equations.

### Divide and Conquer strategy.

The second experiment in `divide_and_conquer_ising1D_anyN.m` continues with this idea of optimizing
the variational parameters of a given wavefunction by finding the set of values that minimize the energy.

In this case, however, instead of finding them analytically, which would be impossible if the number
of parameters were too large, we use a simple exploration scheme called *Divide and Conquer*.

To use this file:
- **Define the problem parameters**:
    - Number of spins `N`.
    - Number of hidden units `h`.
    - Coupling interaction intensity `J`.
    - Magnetic field intensity `B`.
- **Define divide and conquer parameters**:
    - Initial length of square side `length`.
    - *(Optional)* Define the polynomial region for the exploration with `coords`.
    - Depth of the exploration `n_steps`.

This exploration scheme works well for small systems. Nevertheless, the complexity of the energy landscape
as the problem size increases makes this *greedy* algorithm often miss the global minimum and becomes slow
for larger systems (computational cost O($d\cdot n^P$); where $d$ is the depth, $n$ the number of points
defining the exploration region and $P$ the total number of parameters $P = N + h + N\cdot h$).

## Jupyter Notebooks

1. Open the `Jupyter Notebooks` folder.
2. Ensure that you have Jupyter Notebook installed on your system.
3. Launch Jupyter Notebook and navigate to the `Jupyter Notebooks` folder.
4. Run the Jupyter notebooks in the specified order, as mentioned in the experiment documentation.
5. Make sure to provide the required input data or modify the notebook parameters as needed.
6. Save the results and any generated output files in the appropriate location.

This folder contains the notebooks used to generate the figures in the thesis report.
The simulations in these files can be very slow to finish, and should be done multiple times, 
averaging the results to avoid statistical errors in the results.
