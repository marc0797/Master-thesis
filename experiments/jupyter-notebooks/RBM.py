### RBM Class Definition

import numpy as np

class RBM:
    def __init__(self, num_visible, num_hidden, sigma=0.01):
        self.n_visible = num_visible
        self.n_hidden = num_hidden
        self.initialise_parameters(sigma)

    def initialise_parameters(self,sigma):
        ## Initializes parameters as random complex valued gaussians
        b = np.random.randn(self.n_visible) + 1j*np.random.randn(self.n_visible)
        c = np.random.randn(self.n_hidden) + 1j*np.random.randn(self.n_hidden)
        W = (np.random.randn(self.n_hidden,self.n_visible) + 
                1j*np.random.randn(self.n_hidden,self.n_visible))
        self.params = sigma*np.concatenate((b,c,W.ravel()))

    # Define b, c, W as properties in order to avoid accidentally modifying them
    # while getting their values
    @property
    def b(self):
        return self.params[:self.n_visible]

    @property
    def c(self):
        return self.params[self.n_visible:self.n_visible+self.n_hidden]

    @property
    def W(self):
        return np.reshape(self.params[self.n_visible+self.n_hidden:], 
                          (self.n_hidden, self.n_visible))
    def spin_basis(self,N):
        result = []
        for i in range(2 ** N - 1,-1,-1):
            binary_str = format(i, f'0{N}b')
            binary_array = [int(bit) for bit in binary_str]
            result.append(binary_array)
        return result
    
    def log_cosh(self,theta):
        ## Numerically stable log(cosh(theta))
        MAX = 50
        if np.abs(theta) > MAX:
            return np.abs(theta) - np.log(2)
        return np.log(np.cosh(theta))
    
    def wavefunction(self):
        ## Projection of the wavefunction onto spin basis
        log_p = []
        for s in self.spin_basis(self.n_visible):
            theta = self.c + self.W @ s
            exp_b = self.b @ s
            prod_cosh = np.sum([self.log_cosh(theta_i) for theta_i in theta])
            log_p.append(exp_b + prod_cosh)

        # Calculate normalized vector from log_p
        log_max = np.max(log_p)
        scaled_p = np.exp(log_p - log_max)
        scaled_norm = np.sqrt(scaled_p.conj() @ scaled_p)
        return scaled_p/scaled_norm
    
    def p_ratio(self, s1, s2):
        ## Probability ratio between state s2 and reference state s1
        theta1 = self.c + self.W @ s1
        theta2 = self.c + self.W @ s2
        p1 = [self.log_cosh(theta_i) for theta_i in theta1]
        p2 = [self.log_cosh(theta_i) for theta_i in theta2]
        log_diff = self.b @ (s2 - s1) + sum(p2) - sum(p1)
        return np.exp(log_diff)
        
    def p_ratios(self, s1, s2):
        ## Probability ratios between list of states s2 and reference state s1
        return [self.p_ratio(s1, s) for s in s2]