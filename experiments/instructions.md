# Experiment Reproduction Instructions

To reproduce the experiments in this folder, please follow the steps below:

## Python files

The python files in the library contain some useful functions used later in the
notebooks. This is done in order to have "cleaner" notebooks.

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

### Simulated Annealing

With simulated annealing, we should be able to find the optimal set of parameters that minimizes the
energy of the variational wavefunction. This global convergence, however, is only ensured under some
very strict conditions:

1. The starting temperature $T_0$ should be very large compared to the characteristic energy of the system,
and should be decreased infinitessimally at each iteration, for infinitely many iterations.
2. In between temperature modifications, we must let the system thermalize for $N_{\text{thermalization}}\to\infty$.

Since we cannot run this algorithm for an infinite amount of time, in `Ising 1D with Simulated Annnealing.ipynb`
we explore the results of limiting the simulated annealing scheme to a certain number of iterations,
$\gamma$ factors, etc. For some given problem parameters.

### Stochastic Reconfiguration

Stochastic reconfiguration (SR) is a modified version of *gradient descent*, commonly seen in machine learning
problems, that is specifically used for optimizing the parameters in a variational monte carlo problem.

The first notebook that you should see is `Ising 1D with Stochastic Reconfiguration.ipynb`, in which the
basic idea of the method is explained, and a simple example is done for $N$ = 4 spins. 

(*Optional*) There is a statistical analysis of the autocorrelation time of the Markov chain used in the
Monte Carlo sampling scheme of the method. If you want to try it yourself, in Jupyter notebook, set the
corresponding cells to `code` type and run them for different values and problem sizes, try to find what
should be the value of the `block_size` used in SR.

The second and third notebooks (`Ising 1D with SR - N plots.ipynb` and `Ising 1D with SR - B field plots.ipynb`) 
are used to plot the relative error of the method compared to the analytical results
of the problem and the magnetizations of the obtained RBM state, respectively.

**Usage:** Set the path to save the .txt files and the figures. These notebooks should be run multiple times
to obtain results without statistical errors. Ideally, they can also can be written in a python `.py` file
in order to run it in a HPC cluster, with a queue system.
