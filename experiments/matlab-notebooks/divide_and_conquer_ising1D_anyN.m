%% Ising 1D with N spins (BOSONS, h hidden units)
%B
clear
clc

%% Initialization
% Define constant values
N = 4;
h = 1;
J = 2;
B = 1;

% Start with a very simple grid (4 values)
length = 3;
top_left    = -1 + 1j;
top_right   = +1 + 1j;
bot_left    = -1 - 1j;
bot_right   = +1 - 1j;

coords = 0.5*length*ones(N+h+N*h,1).*[top_left,top_right,bot_left,bot_right];
n_steps = 20;

%% Define wavefunction and Hamiltonian
H = buildHamiltonian(N,J,B);

E_min = inf;

energies = eig(H);

fprintf("Exact ground state energy %.3f\n",energies(1))

%% Start divide and conquer strategy

step = 0;

% Initialize parameter vectors and length
params = coords;

tic
while step < n_steps
    % Generate all possible combinations of parameters
    params_combs = my_combvec(params);

    % Iterate through each combination
    for i = 1:size(params_combs,1)
        % Evaluate energy for the combination
        PSI = psi(N,h,params_combs(i,:));
        Norm = PSI'*PSI;
        Energy = PSI'*H*PSI/Norm;

        % Check if current value is smaller than the minimum
        if Energy < E_min
            E_min = Energy;
            min_config = params_combs(i,:);
        end
    end
    % Change vectors according to optimal indexes
    length = length/2;
    params = min_config' + length*coords;

    % Increase step
    step = step + 1;

    % Print results
    if mod(step,5) == 0
        fprintf("\nd = %d\n",step)
        toc
        fprintf("Minimum energy %f\n",E_min)
        fprintf("Relative error with the energy %e\n",...
            (energies(1)-E_min)/energies(1)*100)
    end
end
%% Print results

% fprintf("Minimum parameters %f %+fj\n",real(min_config),...
%     imag(min_config))
% 
% % Optimal psi
% PSI = psi(N,h,min_config);


%% Function to build Hamiltonian

function H = buildHamiltonian(N, J, B)
    % Initialize Pauli matrices
    sigma_x = [0 1; 1 0];
    sigma_z = [1 0; 0 -1];
    
    % Initialize Hamiltonian
    H = zeros(2^N, 2^N);

    % Loop through each term in the Hamiltonian
    for i = 1:N-1
        % Add Ising interaction term
        H = H - J * kron(kron(eye(2^(i-1)), sigma_z),  kron(sigma_z, eye(2^(N-i-1))));
    end
    for i = 1:N
        % Add transverse field term
        H = H - B * kron(kron(eye(2^(i-1)), sigma_x), eye(2^(N-i)));
    end

    % Add periodic boundary conditions (only if N > 2)
    H = H - J * kron(sigma_z, kron(eye(2^(N-2)), sigma_z));
end

%% Wavefunction

function PSI = psi(N,h,params)
    % Generate all spins
    spins = dec2bin(0:2^N-1) - '0';

    % Define RBM parameters (for bosons)
    B = params(1);
    C = params(2:2+h-1);
    W = params(2+h:end);

    exp_term = exp(B*sum(spins,2));
    cosh_term = 1.0;
    for i = 1:h
        cosh_term = 2*cosh_term.*cosh(C(i) + W(i)*sum(spins,2));
    end
    % Return wavevector
    PSI = exp_term.*cosh_term;
    PSI = PSI/sqrt(PSI'*PSI);
end

%% Combvec
function out = my_combvec(params)
    n = size(params,1);
    if n==0
        out = [];
        return
    end
    F = cell(1, n);
    for i = 1:n
        F{i} = params(i,:);
    end
    [F{1:n}] = ndgrid(F{:});
    for i=n:-1:1
        G(:,i) = F{i}(:);
    end
    out = G;
end