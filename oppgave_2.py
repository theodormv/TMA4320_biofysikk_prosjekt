import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from utilities.probs import p_plus, p_minus

def V_1(x, k):
    return k

def V_2(x, k):
    return -k * x

def random_walk(V):
    k_b = sp.constants.Boltzmann
    n_particles = 100
    n_steps = 100
    h = 1
    dt = 1
    T = 273.15 + 37
    beta = (T * k_b)**-1
    k = 1/beta    
    dist_list = list()

    for i in range(n_steps):
        particles = np.arange(1, n_particles+1, 1)
        uniform_dist = np.random.uniform(0, 1, n_particles)
        x_0 = particles
        go_left = uniform_dist <= p_minus(x_0, beta, k, V)
        go_right = uniform_dist >= 1-p_plus(x_0, beta, k, V)
        stay = (~go_left & ~go_right)
        particles[go_left] -= 1
        particles[go_right] += 1
        dist_list.append(particles)

        if i % 10 == 0:
            print(f'p_plus: {np.average(p_plus(x_0, beta, k, V))} | p_minus: {np.average(p_minus(x_0, beta, k, V))}')
    
    return np.array(dist_list).flatten()

'''
oppg_1_a_verdier = random_walk(V_1)
plt.hist(oppg_1_a_verdier, density=True)
plt.show()
'''

oppg_2_a_verdier = random_walk(V_2)
plt.hist(oppg_2_a_verdier, density=True)
plt.show()