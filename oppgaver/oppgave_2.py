import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from utilities.probs import p_plus, p_minus

def V_1(x):
    return k

def V_2(x):
    return - k * x


class random_walk():
    def __init__(self, beta_k_ratio : int, V):
        k_b = sp.constants.Boltzmann
        self.n_particles = 100
        self.n_walks = 200
        self.V = V
        T = 273.15 + 37
        self.beta = (T * k_b)**-1
        self.k = beta_k_ratio * (T * k_b)
        global k
        k = self.k
        self.beta_k_ratio = beta_k_ratio
        self.oppg_a_vals = list()
        self.oppg_b_vals = list()

    def gen_walk_masks(self, particles):
        uniform_dist = np.random.uniform(0, 1, self.n_particles)
        x_0 = particles
        go_left = uniform_dist <= p_minus(x_0, self.beta, self.V)
        go_right = uniform_dist >= 1-p_plus(x_0, self.beta, self.V)
        return go_left, go_right

    def oppg_a_step(self):
        particles = np.arange(0, self.n_particles, 1)
        for j in range(self.n_walks):
            go_left, go_right = self.gen_walk_masks(particles)
            particles[go_left] -= 1
            particles[go_right] += 1
        self.oppg_a_vals.append(particles)
    
    def oppg_b_step(self):
        particles = np.arange(0, self.n_particles, 1)
        for j in range(self.n_walks):
            go_left, go_right = self.gen_walk_masks(particles)
            movement_order = np.random.permutation(np.arange(0,self.n_particles))
            for idx in movement_order:
                if go_left[idx] == True:
                    if np.any(particles[idx] -1 == particles):
                        continue
                    else:
                        particles[idx] -= 1
                if go_right[idx] == True:
                    if np.any(particles[idx] + 1 == particles):
                        continue
                    else:
                        particles[idx] += 1
        self.oppg_b_vals.append(particles)


n_steps = 100

walker_V1 = random_walk(beta_k_ratio=1000, V=V_1)
walker_V2 = random_walk(beta_k_ratio=1000, V=V_2)

for j in range(n_steps):
    walker_V1.oppg_a_step()
    walker_V1.oppg_b_step()
    walker_V2.oppg_a_step()
    walker_V2.oppg_b_step()
    if j % 10 == 0:
        print(f'Iteration {j} of {n_steps}')


for walker, potential_info in {walker_V1:('V1', 'k'), walker_V2: ('V2', '-kx')}.items():
    potential_type, potential_eq = potential_info


    verdier = np.array(walker.oppg_a_vals).flatten()
    plt.hist(verdier, density=True, bins=walker.n_particles//3)
    plt.title(f'Partikkelfordeling uten interaksjoner og $V(x) = {potential_eq}$')
    plt.savefig(f'figures\\2a_{potential_type}_{walker.beta_k_ratio}ratio.jpg')

    plt.clf()
    

    verdier = np.array(walker.oppg_b_vals).flatten()
    plt.hist(verdier, density=True, bins=walker.n_particles//3)
    plt.title(f'Partikkelfordeling med interaksjoner og $V(x) = {potential_eq}$')
    plt.savefig(f'figures\\2b_{potential_type}_{walker.beta_k_ratio}ratio.jpg')


    plt.clf()
