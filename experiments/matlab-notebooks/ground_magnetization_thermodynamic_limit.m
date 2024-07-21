%% Ising 1D Exact Magnetization of the Ground State:
% Using the ground state found with Jordan-Wigner and Bogoliubov
% transformations, we want to calculate its exact magnetization:
%
%   m_z = \langle\sigma^z\rangle
%
%   m_x = \langle\sigma^x\rangle
%
% which can be calculated as the mean of the products:
%
%   m_z = \frac{1}{L}\sum_{i=1}^{L}\ev*{\hat{\sigma}^z}{\emptyset_\gamma}
%
%   m_x = \frac{1}{L}\sum_{i=1}^{L}\ev*{\hat{\sigma}^x}{\emptyset_\gamma}
%
% Also, I have found an approximation for large L of the previous
% expression, that simplifies them to an integral:
%
%   m_z = (1 - \frac{B^2}{J^2})^{1/8}/2 
%
%   m_x = \frac{1}{2\pi}\int_0^{\pi}(\frac{B + Jcos(k)}{\sqrt{J^2 +...
%           B^2 + 2JBcos(k)}} dk
%
% So, in order to verify the validity of this approximation, we will first
% plot the error as a function of L, and then plot the magnetization as a
% function of B/J to observe the phase transition at B = J.

clear
close all
clc

%========================== HARD LIMIT: N < 20 ==========================%

%% Error plot as a function of N
%{
N_vec = 2:2:12;
mz_errors = 0*N_vec;
mx_errors = 0*N_vec;

for i = 1:length(N_vec)
    % Find Hamiltonian and exact ground state
    H = buildHamiltonian(N_vec(i),J,B);
    
    [state,energy] = eigs(H,1,'smallestreal');
    
    % Find exact mean magnetization, both longitudinal and transversal
    [Mz, Mx] = buildMagnetization(N_vec(i));
    
    exact_mz = sqrt(state'*Mz^2*state);
    exact_mx = sqrt(state'*Mx^2*state);

    % calculate error
    rel_error_mz = abs((exact_mz - aprox_mz)/exact_mz);
    rel_error_mx = abs((exact_mx - aprox_mx)/exact_mx);

    mz_errors(i) = rel_error_mz;
    mx_errors(i) = rel_error_mx;
end

% Plot
figure(1);
loglog(2:2:12,mz_errors,'*b')
hold on
loglog(2:2:12,mx_errors,'*r')
title("Relative error of the approximation")
xlabel("L")
ylabel("error")
legend("m_z","m_x")
%}

%% Magnetization as a function of B/J
tic
% fix J and N (suppose N -> inf)
J = 1;

x_vec = linspace(-1,1,20);
z_vec = (2-0.25)*x_vec.^3 + 0.25*x_vec;
B_vec = 10.^z_vec;

L = 5;
exmz_vec = zeros(L,length(B_vec));
exmx_vec = zeros(L,length(B_vec));
exchi_vec = zeros(L,length(B_vec));
exchiz_vec = zeros(L,length(B_vec));
mz_vec = 0*B_vec;
mx_vec = 0*B_vec;
chi_vec = 0*B_vec;

for i = 1:length(B_vec)
    % Calculate magnetization integrals
    [aprox_mz, aprox_mx] = magnetizationIntegrals(J,B_vec(i));
    
    mz_vec(i) = aprox_mz;
    mx_vec(i) = aprox_mx;

    for j = 1:L
        % Find exact magnetization matrices
        [Mz, Mx] = buildMagnetization(4*j);

        % Find Hamiltonian and exact ground state
        H = buildHamiltonian(4*j,J,B_vec(i));
    
        [state,energy] = eigs(H,1,'smallestreal');

        exmz_vec(j,i) = sqrt(state'*Mz^2*state);
        exmx_vec(j,i) = sqrt(state'*Mx^2*state);
        exchi_vec(j,i) = 4*j*(state'*Mx^2*state - (state'*abs(Mx)*state)^2);
        exchiz_vec(j,i) = 4*j*(state'*Mz^2*state - (state'*abs(Mz)*state)^2);
    end
end

%%

figure(2);
for k = 1:5
    semilogx(B_vec,exmz_vec(k,:),'-o','LineWidth',2.5,...
        'MarkerSize',6,'DisplayName',['N = ',num2str(4*k)])
    hold on
end
semilogx(B_vec,mz_vec,'-k','LineWidth',2.5,'DisplayName','N = \infty')
title("\langle m_z\rangle as a function of B/J")
xlabel("B/J")
ylabel("m_z")
legend('show')

figure(3);
for k = 1:5
    semilogx(B_vec,exmx_vec(k,:),'-o','LineWidth',2.5,...
        'MarkerSize',6,'DisplayName',['N = ',num2str(4*k)])
    hold on
end
semilogx(B_vec,mx_vec,'-k','LineWidth',2.5,'DisplayName','N = \infty')
title("\langle m_x\rangle as a function of B/J")
xlabel("B/J")
ylabel("m_x")
legend('show')

figure(4);
for k = 1:5
    semilogx(B_vec,exchi_vec(k,:),'-o','LineWidth',2.5,...
        'MarkerSize',6,'DisplayName',['N = ',num2str(4*k)])
    hold on
end
semilogx(B_vec,chi_vec,'-k','LineWidth',2.5,'DisplayName','N = \infty')
title("Magnetic (\chi_x) susceptibility as a function of B/J")
xlabel("B/J")
ylabel("\chi_x")
legend('show')

figure(5);
for k = 1:5
    semilogx(B_vec,exchiz_vec(k,:),'-o','LineWidth',2.5,...
        'MarkerSize',6,'DisplayName',['N = ',num2str(4*k)])
    hold on
end
semilogx(B_vec,chi_vec,'-k','LineWidth',2.5,'DisplayName','N = \infty')
title("Magnetic (\chi_z) susceptibility as a function of B/J")
xlabel("B/J")
ylabel("\chi_z")
legend('show')

toc
%% Functions

% A function to find the ground state of the Ising 1D chain
function H = buildHamiltonian(N, J, B)
    % Initialize Pauli matrices
    sigma_x = sparse([0 1; 1 0]);
    sigma_z = sparse([1 0; 0 -1]);
    
    % Initialize Hamiltonian
    H = - J * kron(kron(speye(1), sigma_z),  kron(sigma_z, speye(2^(N-2))));

    % Loop through each term in the Hamiltonian
    for i = 2:N-1
        % Add Ising interaction term
        H = H - J * kron(kron(speye(2^(i-1)), sigma_z),...
            kron(sigma_z, speye(2^(N-i-1))));
    end
    for i = 1:N
        % Add transverse field term
        H = H - B * kron(kron(speye(2^(i-1)), sigma_x), speye(2^(N-i)));
    end

    % Add periodic boundary conditions (only if N > 2)
    if N >= 2
        H = H - J * kron(sigma_z, kron(speye(2^(N-2)), sigma_z));
    end
end

% A function to find the mean magnetization
function [Mz, Mx] = buildMagnetization(N)
    % Initialize Pauli matrices (only the nonzero terms)
    sigma_x = sparse([0 1; 1 0]);
    sigma_z = sparse([1 0; 0 -1]);

    Mz = kron(kron(speye(1), sigma_z), speye(2^(N-1)));
    Mx = kron(kron(speye(1), sigma_x), speye(2^(N-1)));
    % Loop through each lattice site
    for i = 2:N
        % Add longitudinal term
        Mz = Mz + kron(kron(speye(2^(i-1)), sigma_z), speye(2^(N-i)));
        
        % Add transversal term
        Mx = Mx + kron(kron(speye(2^(i-1)), sigma_x), speye(2^(N-i)));
    end
    Mz = Mz/N;
    Mx = Mx/N;
end

% A function to calculate the magnetization integrals
function [Mz, Mx] = magnetizationIntegrals(J,B)
    % Define integrals
    f = @(k) (B + J*cos(k))./sqrt(J^2 + B^2 + 2*J*B*cos(k));

    % Calculate magnetizations
    if B > J
        Mz = 0;
    else
        Mz = (1 - B^2/J^2)^(1/8);
    end
    Mx = integral(f,0,pi)/pi;
end