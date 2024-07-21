%% Ising 1D Exact Ground State energy:
% k-basis for discrete Fourier transformation of the fermionic operators:
% If the number of lattice sites L is even, we can define the positive
% k-value basis as:
%
%   K^+_{p=1} = \{k = 2n\pi/L, with n = 1, ..., L/2 - 1\}
%   K^+_{p=0} = \{k = (2n - 1)\pi/L, with n = 1, ..., L/2\}
%
% The ground state energy for the L even, p = 0 (ABC) case is 
% 
%   E^{ABC}_0 = -\sum_{k\in K^+_{p=0}}\epsilon_k
%
% In the L even, p = 1 (PBC) case, we need to add an extra term
%
%   E^{PBC}_0 = - 2J - \sum_{k\in K^+_{p=1}}\epsilon_k
%
% However, in the case where the number of lattice sites is odd, we have to
% keep the whole expression for K_{p=0,1}:
%
%   K'_{p=1} = \{k = 2n\pi/L, with n = -(L-1)/2, ..., 0, ..., (L-1)/2\}
%   K'_{p=0} = \{k = (2n+1)\pi/L, with n = -(L-1)/2, ..., (L-1)/2\}
%
% p being the parity of the number of fermions (N), p = 0 even, p = 1 odd.
% Notice that in the K_{p=1}, we can redefine the values of k to form BCS
% pairs as in the previous L even case:
%
%   K'^+_{p=1} = \{k = 2n\pi/L, with n = 1, ..., (L-1)/2\}
%
% while the k = 0 case contributes -|J - B| to the final energy
%
%   E^{PBC}_0 = - |J - B| - \sum_{k\in K'^+_{p=1}}\epsilon_k
% 
% The k values in the p = 0 case can also be written to form BCS pairs 
% as:
%
%   K'^+_{p=0} = \{k = (2n + 1)\pi/L, with n = 0, ..., (L-3)/2\}
% 
% In this case, however, since we have an even number of fermions, the
% contribution of n_\pi should be 0.
% 
%   E^{ABC}_0 = -\sum_{k\in K'^+_{p=0}}\epsilon_k
%
% We will see in this code, that the true ground state energy for the L
% even case corresponds to antiperiodic boundary conditions or p = 0, while
% in the L odd case it is the opposite.
% If we assume N = L, then their parity matches and we always end up with
% the true ground state.

clear
close all
clc

%% Thermodynamic limit, L even case

L = 10000;  % L -> inf
J = 1;      % We fix the interaction term, and vary the magnetic term

delta_E = zeros(5,1);
for B = 0:4
    delta_E(B+1) = ising_ground_energy(L,"odd",J,B) ...
        - ising_ground_energy(L,"even",J,B);
end

% Plot of the difference between periodic and antiperiodic boundary states
figure(1)
plot(0:4,delta_E,'-b')
title("Energy difference between PBC and ABC, L even")
xlabel("B/J")
ylabel("$E^{PBC}_0 - E^{ABC}_0$",'Interpreter','latex')
text(3,2,"2(B-J)",...
    "HorizontalAlignment","center","FontSize",12)
axis([0 4 -0.5 6.5])
xlims = get(gca,'xlim');
ylims = get(gca,'ylim');
annotation("arrow",[2.875 2.875]/(xlims(2)-xlims(1)),[2.9 4.2]/(ylims(2)-ylims(1)))

% Plot of the ferromagnetic phase
L_vec = 4:4:40;
i = 1;
ferr_delta_E = zeros(length(L_vec),1);
for L_i = L_vec
    ferr_delta_E(i) = ising_ground_energy(L_i,"odd",J,0.5) ...
        - ising_ground_energy(L_i,"even",J,0.5);
    i = i + 1;
end
figure(2)
semilogy(L_vec,ferr_delta_E,'-ob','MarkerFaceColor','blue','MarkerSize',4)
title("Ferromagnetic phase: B = J/2, L even")
xlabel("L")
ylabel("$E^{PBC}_0 - E^{ABC}_0$",'Interpreter','latex')
axis([0 45 1e-16 1])

% Plot of the critical point
L_vec = [10:10:100 1000 10000];
i = 1;
crit_delta_E = zeros(length(L_vec),1);
for L_i = L_vec
    crit_delta_E(i) = ising_ground_energy(L_i,"odd",J,J) ...
        - ising_ground_energy(L_i,"even",J,J);
    i = i + 1;
end
figure(3)
loglog(L_vec,crit_delta_E,'ob','MarkerFaceColor','blue','MarkerSize',4)
hold on
loglog(L_vec,pi./(2*L_vec),'-r')
title("Critical point: B = J, L even")
xlabel("L")
ylabel("$E^{PBC}_0 - E^{ABC}_0$",'Interpreter','latex')
text(3e2,2e-2,["Critical point:","\pi/(2L)","\downarrow"],...
    "HorizontalAlignment","center","FontSize",12)
axis([0 1.5e4 1e-4 1])

%% Thermodynamic limit, L odd case

L = 10001;  % L -> inf
J = 1;      % We fix the interaction term, and vary the magnetic term

delta_E = zeros(5,1);
for B = 0:4
    delta_E(B+1) = ising_ground_energy(L,"even",J,B) ...
        - ising_ground_energy(L,"odd",J,B);
end

% Plot of the difference between antiperiodic and periodic boundary states
figure(4)
plot(0:4,delta_E,'-b')
title("Energy difference between ABC and PBC, L odd")
xlabel("B/J")
ylabel("$E^{ABC}_0 - E^{PBC}_0$",'Interpreter','latex')
text(3,2,"2(B-J)",...
    "HorizontalAlignment","center","FontSize",12)
axis([0 4 -0.5 6.5])
xlims = get(gca,'xlim');
ylims = get(gca,'ylim');
annotation("arrow",[2.875 2.875]/(xlims(2)-xlims(1)),[2.9 4.2]/(ylims(2)-ylims(1)))

% Plot of the ferromagnetic phase
L_vec = 3:4:39;
i = 1;
ferr_delta_E = zeros(length(L_vec),1);
for L_i = L_vec
    ferr_delta_E(i) = ising_ground_energy(L_i,"even",J,0.5) ...
        - ising_ground_energy(L_i,"odd",J,0.5);
    i = i + 1;
end
figure(5)
semilogy(L_vec,ferr_delta_E,'-ob','MarkerFaceColor','blue','MarkerSize',4)
title("Ferromagnetic phase: B = J/2, L odd")
xlabel("L")
ylabel("$E^{ABC}_0 - E^{PBC}_0$",'Interpreter','latex')
axis([0 45 1e-16 1])

% Plot of the critical point
L_vec = [11:10:101 1001 10001];
i = 1;
crit_delta_E = zeros(length(L_vec),1);
for L_i = L_vec
    crit_delta_E(i) = ising_ground_energy(L_i,"even",J,J) ...
        - ising_ground_energy(L_i,"odd",J,J);
    i = i + 1;
end
figure(6)
loglog(L_vec,crit_delta_E,'ob','MarkerFaceColor','blue','MarkerSize',4)
hold on
loglog(L_vec,pi./(2*L_vec),'-r')
title("Critical point: B = J, L odd")
xlabel("L")
ylabel("$E^{ABC}_0 - E^{PBC}_0$",'Interpreter','latex')
text(3e2,2e-2,["Critical point:","\pi/(2L)","\downarrow"],...
    "HorizontalAlignment","center","FontSize",12)
axis([0 1.5e4 1e-4 1])
%% Functions

function K = k_values(L,N)
    % Function that calculates the values that should be used in the
    % Fourier transformation of the fermionic operators.

    if mod(L,2) == 0 && N == "even"
        % If L is even, and N is even, calculate K^+_{p = 0}
        K = (2*(1:L/2) - 1)*pi/L;
    elseif mod(L,2) == 0 && N == "odd"
        % If L is even, and N is odd, calculate K^+_{p = 1} instead
        K = 2*(1:L/2-1)*pi/L;
    elseif mod(L,2) == 1 && N == "even"
        % If L is odd, and N is even, calculate K_{p = 0}
        K = (2*(0:(L-3)/2) + 1)*pi/L;
    else
        % If L is odd, and N is odd, calculate K_{p = 1}
        K = 2*(1:(L-1)/2)*pi/L;
    end
end

function E = ising_ground_energy(L,N,J,B)
    % Function that calculates the exact ground state energy of the
    % quantum uniform Ising model

    % Find the k-values
    K = k_values(L,N);

    % The ground state energy will be the sum of all \epsilon_{k-}
    E = -2*J*sum(sqrt((cos(K) - B/J).^2 + sin(K).^2));
    if mod(L,2) == 0 && N == "odd"
        % If L is even and N is odd, then we need to add an extra term
        E = E - 2*J;
    elseif mod(L,2) == 1 && N == "odd"
        % If L is odd, we need to divide the energy contribution in half
        E = E - abs(J - B) - (J - B);
    elseif mod(L,2) == 1 && N == "even"
        % If L is odd, and N is even, we also need to add an extra term
        E = E;
    end
end