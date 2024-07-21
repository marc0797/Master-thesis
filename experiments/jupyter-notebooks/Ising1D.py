## Functions for exact diagonalization and magnetizations

import numpy as np
from scipy.sparse import kron, csr_matrix
from scipy.sparse.linalg import eigsh

def buildHamiltonian_sparse(N, J, B):
    # Sparse Hamiltonian for the Transverse Field Ising model
    
    # Define the Pauli matrices
    Id = csr_matrix([[1, 0], [0, 1]], dtype=float)
    sigma_z = csr_matrix([[1, 0], [0, -1]], dtype=float)
    sigma_x = csr_matrix([[0, 1], [1, 0]], dtype=float)

    # Coupling term
    first_term = [Id] * N
    first_term[0] = sigma_z
    first_term[1] = sigma_z

    # Magnetic term
    second_term = [Id] * N
    second_term[0] = sigma_x

    # Initialize sparse zero matrix
    H = csr_matrix((2**N, 2**N), dtype=float)

    # Build Hamiltonian
    for i in range(N):
        # Compute the Kronecker product for the current terms
        first_term_kronecker = first_term[0]
        for mat in first_term[1:]:
            first_term_kronecker = kron(first_term_kronecker, mat, format='csr')
        
        second_term_kronecker = second_term[0]
        for mat in second_term[1:]:
            second_term_kronecker = kron(second_term_kronecker, mat, format='csr')
        
        # Add to Hamiltonian
        H += - J*first_term_kronecker - B*second_term_kronecker
        
        # Circular shift
        first_term = [first_term[-1]] + first_term[:-1]
        second_term = [second_term[-1]] + second_term[:-1]
    
    return H

def diagonalizeHamiltonian_sparse(H,T=1):
    # Diagonalize the Hamiltonian to find the ground state.
    # Use eigsh to find the smallest eigenvalue and the corresponding eigenvector
    eigenvalues, eigenvectors = eigsh(H, k=T, which='SA')

    energies = eigenvalues[:T]
    states = eigenvectors[:,:T]
    
    return energies, states

def groundEnergy(N,J,B):
    # Excitation energy of the Hamiltonian given J and B
    N_half = int(N/2)
    # Define k vectors
    k_vec = [np.pi*(2*n - 1)/N for n in range(1,N_half+1)]

    cos_term = (np.cos(k_vec) - B/J)**2
    sin_term = np.sin(k_vec)**2
    return -2*J*np.sum(np.sqrt(cos_term + sin_term))

def buildMagnetizations_sparse(N):
    # Sparse Magnetizations for the Transverse Field Ising model
    
    # Define the Pauli matrices
    Id = csr_matrix([[1, 0], [0, 1]], dtype=float)
    sigma_z = csr_matrix([[1, 0], [0, -1]], dtype=float)
    sigma_x = csr_matrix([[0, 1], [1, 0]], dtype=float)

    # Coupling term
    mz_term = [Id] * N
    mz_term[0] = sigma_z

    # Magnetic term
    mx_term = [Id] * N
    mx_term[0] = sigma_x

    # Initialize sparse zero matrix
    Mz = csr_matrix((2**N, 2**N), dtype=float)
    Mx = csr_matrix((2**N, 2**N), dtype=float)

    # Build Magnetizations
    for i in range(N):
        # Compute the Kronecker product for the current terms
        mz_term_kronecker = mz_term[0]
        for mat in mz_term[1:]:
            mz_term_kronecker = kron(mz_term_kronecker, mat, format='csr')
        
        mx_term_kronecker = mx_term[0]
        for mat in mx_term[1:]:
            mx_term_kronecker = kron(mx_term_kronecker, mat, format='csr')
        
        # Add to Magnetization matrices
        Mz += mz_term_kronecker/N
        Mx += mx_term_kronecker/N
        
        # Circular shift
        mz_term = [mz_term[-1]] + mz_term[:-1]
        mx_term = [mx_term[-1]] + mx_term[:-1]
    
    return Mz, Mx

def FubiniStudy(psi,phi):
    return np.arccos(np.sqrt(psi.conj().T @ phi * phi.conj().T @ psi))