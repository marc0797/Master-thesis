%% Ising 1D with N = 2 spins
%B
clear
clc

%% Initialization
% Define variables symbolically
syms s1 s2 b J B

% assumptions
assume(b,"Real")
J = 2;
B = 1;

%% Define wavefunction and Hamiltonian
H = [   -2*J   -B   -B    0;
          -B  2*J    0   -B;
          -B    0  2*J   -B;
           0   -B   -B -2*J];
psi(s1,s2) = exp(b*s1 + b*s2);

PSI = [psi(1,1); psi(1,0); psi(0,1); psi(0,0)];

%% Solve for the minimum mean energy

% Evaluate mean energy
N = conj(PSI')*PSI;
E = conj(PSI')*H*PSI/N;

% Differentiate with respect to the weights
dE = diff(E,b) == 0;

% Solve system of equations
solutions = solve(dE)


%% Calculate minimum energy

Emin = inf;
for i = 1:length(solutions)
    %psi_min(s1,s2) = (solutions.a1(i) + 1i*solutions.b1(i))^s1*...
    %    (solutions.a2(i) + 1i*solutions.b2(i))^s2;
    psi_min(s1,s2) = exp(solutions(i)*(s1+s2));

    PSI_min = [psi_min(1,1), psi_min(1,0), psi_min(0,1), psi_min(0,0)];

    N_min = conj(PSI_min)*PSI_min';
    E_min = conj(PSI_min)*H*PSI_min'/N_min;
    eval(E_min)
    if eval(E_min) < Emin
        Emin = eval(E_min);
        b1_min = solutions(i);
    end
end
energies = eig(H);

fprintf("Exact ground state energy %f\n",energies(1))
fprintf("Minimum energy %.2f\n",Emin)
fprintf("Relative error with the energy %f%%\n",...
    (energies(1)-Emin)/energies(1)*100)
fprintf("Minimum parameters b1 = b2 = b = %f %+fj\n",real(b1_min),...
    imag(b1_min))