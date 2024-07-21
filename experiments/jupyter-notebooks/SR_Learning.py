### Useful functions for Stochastic Reconfiguration learning

import numpy as np
from copy import deepcopy
from Ising1D import buildMagnetizations_sparse

def local_energy(s, psi, pbc=True,J=2,B=1):
    ## Local energy of TFI model
    # Interaction term
    couplings = (s[:-1]==s[1:])*2-1
    e_interaction = -J*sum(couplings)

    # Add Periodic Boundary Conditions
    if pbc:
        e_interaction -= J*((s[-1]==s[0])*2-1)
    
    # Transverse field
    states_with_flip = [flip(s, i) for i in range(len(s))]
    e_field = -B*sum(psi.p_ratios(s, states_with_flip))

    return e_interaction + e_field 
    
def flip(s, i):
    ## Flips i-th spin of state s
    s_flip = deepcopy(s)
    s_flip[i] = 1 - s_flip[i]
    return s_flip

def bin_averages(x,bs):
    ## Bins time-series x into bs chunks and takes means
    nb = len(x)//bs
    bin_avg = [np.mean(x[block*bs:(block+1)*bs]) for block in range(nb)]
    return np.array(bin_avg)

def variational_derivative(s, psi):
    ## Variational log-derivatives for SR
    theta = psi.c + psi.W @ s
    Ob = s
    Oc = np.tanh(theta)
    Ow = Oc[:, None] @ s[None, :]
    return np.concatenate((Ob, Oc, Ow.ravel()))

def covariance(s1, s2):
    ## Computes the covariance between states s1 and s2.
    samples = s1.shape[1]
    m1 = np.mean(s1, axis=1)
    m2 = np.mean(s2, axis=1) if len(s2.shape)>1 else np.mean(s2)
    return (s1 @ s2.T)/samples - m1[:,None] @ m2[None,:]

def sample_block(psi, bs, s0=None, n_flips=1):
    ## Sample bs states according to psi.
    state = np.random.randint(0, 2, psi.n_visible) if s0 is None else s0
    states = []
    for _ in range(bs):
        spin_idx = np.random.randint(0, psi.n_visible, n_flips)
        new_state = flip(state, spin_idx)
        if np.random.random() <= np.abs(psi.p_ratio(state, new_state))**2:
            state = deepcopy(new_state)   # Accept new state   
        states.append(state)
    return states

def stochastic_reconfiguration(psi, J, B, iters=375, lr=1e-2, l2r=1e-3, n_blocks=150, bs=10, n_flips=1, bins=100):
    ## Stochastic reconfiguration function
    thermalise = int(0.1*n_blocks)
    nb = n_blocks - thermalise

    n = psi.n_visible
    Mz, Mx = buildMagnetizations_sparse(n)

    # Empty energies array
    energies = []
    
    # Main loop
    for it in range(iters):
        EL, O = np.zeros(nb, dtype=complex), np.zeros((len(psi.params), nb), dtype=complex)
        states = sample_block(psi, thermalise*bs, n_flips=n_flips)
        state = states[-1]
        for k in range(nb):
            batch = sample_block(psi, bs, s0=state, n_flips=n_flips)
            states += batch
            state = batch[-1]
            EL[k] = local_energy(state, psi,J=J,B=B)
            O[:, k] = variational_derivative(state, psi)
        
        energies.append(EL.mean())
        F = covariance(O.conj(), EL[None,:])   # Gradient
        S = covariance(O.conj(), O)            # Fisher info
        Sinv = np.linalg.pinv(S, rcond=1e-5)   # (pseudo)Inversion
        d_params = lr*Sinv @ F
        psi.params -= d_params.squeeze() + l2r*lr*psi.params

    # RBM ground state energy
    RBM_energy = np.mean(np.real(energies[-bins:]))

    # Calculate the projections of our wavefunction onto the spin basis
    ket_psi = psi.wavefunction()

    RBM_mz = np.sqrt(ket_psi.conj().T @ Mz**2 @ ket_psi)
    RBM_mx = np.sqrt(ket_psi.conj().T @ Mx**2 @ ket_psi)
    RBM_chix = n*(ket_psi.conj().T @ Mx**2 @ ket_psi - (ket_psi.conj().T @ Mx @ ket_psi)**2)
    RBM_chiz = n*(ket_psi.conj().T @ Mz**2 @ ket_psi - (ket_psi.conj().T @ Mz @ ket_psi)**2)

    return RBM_energy, RBM_mz, RBM_mx, RBM_chix, RBM_chiz