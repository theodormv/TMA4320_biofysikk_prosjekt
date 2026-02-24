import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

from utilities.probs import p_minus, p_plus



def V_1(x : np.array):
    k = cfg['alpha']
    N_x = cfg['N_x']
    alpha = cfg['alpha']
    
    #Mapper x til periodisiteten
    L = -(1 - alpha) * N_x
    x = L + (x - L) % N_x
   
    #delt funksjonsuttrykk
    x_less = (0 <= x) & (x <= alpha*N_x)
    x_more = (-(1-alpha)*N_x < x) & (x <= 0)
    x[x_less] = k * x[x_less]/ (alpha*N_x)
    x[x_more] = -k*x[x_more]/((1-alpha)*N_x)
    
    return x

def V_2(x):
    return k


cfg = {
    'alpha': 0.2,
    'T': 273.15 + 37,
    'N_walks': 100,
    'N_x': 100,
    'T_p': 300,
    'N_cycles': 5,
    'N_s': 10,
    'h': 1,
    'beta_k_ratio': 1000,
    'N_x/N_p': 1,
    'b': 2, # Partikkelstørrelse
    'potentials': [V_1, V_2]
}



class ratchet_interaction_walker():
    def __init__(self, cfg : dict):
        k_b = sp.constants.Boltzmann
        self.n_walks = cfg['N_walks']
        self.potentials = cfg['potentials']
        self.T = cfg['T']
        self.beta = (self.T * k_b)**-1
        self.beta_k_ratio = cfg['beta_k_ratio']
        self.k = self.beta_k_ratio * (self.T * k_b)
        global k
        k = self.k
        self.h = cfg['h']
        self.b = cfg['b']
        self.T_p = cfg['T_p']
        self.N_cycles = cfg['N_cycles']
        self.alpha = cfg['alpha']
        self.N_x = cfg['N_x']
        self.N_points = self.N_x * cfg['N_s']
        self.N_particles = self.N_x // cfg['N_x/N_p']

    def gen_walk_masks(self, x_0, V):
        uniform_dist = np.random.uniform(0, 1, self.N_particles)
        go_left = uniform_dist <= p_minus(x_0, self.beta, V)
        go_right = uniform_dist >= 1-p_plus(x_0, self.beta, V)
        return go_left, go_right

    
    def interaction_simulator(self):
        '''Simulerer "Random walk in a ratchet potential with interactions". Returnerer tidssteg og x-posisjon til tilfeldige partikler gjennom simulering'''
        potential_switch_count = -1 # Satt til -1 slik at første 40 steps er ved potensial V_1
        particles = np.arange(0, self.N_points, self.N_points//self.N_particles)

        #Velger tilfeldige partikler som plottes
        plt_particle_idxs = np.random.choice(np.arange(0, self.N_particles-1), size = 5, replace=False)
        particle_movements = list()

        for j in range(self.T_p * 2 * self.N_cycles):

            #Itererer periodisk gjennom potensialer etter T_p tidssteg
            if j % self.T_p == 0:
                potential_switch_count += 1
                V = self.potentials[potential_switch_count % len(self.potentials)]
                print(f'Cycle {j}: Bytter til potensialet på idx {potential_switch_count % len(self.potentials)}')

            go_left, go_right = self.gen_walk_masks(particles, V)
            movement_order = np.random.permutation(np.arange(self.N_particles))
            for curr_particle_idx in movement_order:
                distance_from_particle = np.copy(particles) - particles[curr_particle_idx]
                if go_left[curr_particle_idx]:
                    if np.any((distance_from_particle >= -self.b) & (distance_from_particle < 0)):
                        continue
                    else:
                        #Denne if-statementen sjekker om periodisk randbetingelse er brutt
                        if curr_particle_idx < np.ceil(self.b):
                            boundary = self.N_points - 1 - (self.b - curr_particle_idx) # -1 pga 0-indeksering
                            if np.any(distance_from_particle >= boundary):
                                continue
                    #Denne koden kjører kun dersom ingen if-statements ovenfor er gyldig
                    particles[curr_particle_idx] -= 1

                elif go_right[curr_particle_idx]:
                    if np.any((distance_from_particle <= self.b) & (distance_from_particle > 0)):
                        continue
                    if curr_particle_idx > self.N_points - 1 - self.b:
                        boundary = self.b - (self.N_points - 1 - curr_particle_idx)
                        if np.any(np.abs(distance_from_particle) >= self.N_points - 1 - boundary):
                            continue

                    particles[curr_particle_idx] += 1

                #Syklisk randbetingelse for x
                if particles[curr_particle_idx] == self.N_points:
                    particles[curr_particle_idx] = 0
                
                elif particles[curr_particle_idx] == -1:
                    particles[curr_particle_idx] = self.N_points - 1
            
            particle_movements.append(particles[plt_particle_idxs])

        particle_movements = np.array(particle_movements).T
        return (np.arange(self.N_cycles * self.T_p * 2), particle_movements)
    

if __name__ == '__main__':
    #Bare for å se hvordan potensialet ser ut. Fjern senere.
    walker = ratchet_interaction_walker(cfg)
    t = np.linspace(0, walker.N_points, 2000)
    x = V_1(t)
    plt.plot(t, x)
    plt.show()