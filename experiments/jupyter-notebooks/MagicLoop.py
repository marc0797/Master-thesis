## This file is just for the definition of a parallelized loop

from RBM import RBM
import SR_Learning as sr

def magic_function(n,B):
    # Initialising RBM
    n_visible, n_hidden = n, 2*n
    psi = RBM(n_visible, n_hidden, sigma=0.1)

    RBM_energy, RBM_mz, RBM_mx, RBM_chix, RBM_chiz = sr.stochastic_reconfiguration(psi, 1, B)
    return RBM_energy, RBM_mz, RBM_mx, RBM_chix, RBM_chiz