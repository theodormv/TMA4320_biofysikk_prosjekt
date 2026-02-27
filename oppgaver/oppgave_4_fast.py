import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from numba import njit

from utilities.probs_fast import p_minus, p_plus
np.seterr(over='ignore')

@njit
def V_1(x : np.array):
    
    #Mapper x til periodisiteten
    L = -(1 - alpha) * N_x
    x = L + (x - L) % N_x

    #delt funksjonsuttrykk
    x_less = (0 <= x) & (x <= alpha*N_x)
    x_more = (-(1-alpha)*N_x < x) & (x <= 0)
    x[x_less] = k * x[x_less]/ (alpha*N_x)
    x[x_more] = -k*x[x_more]/((1-alpha)*N_x)
    
    return x

@njit
def V_2(x):
    return np.full_like(x, k, dtype=np.float64)

@njit
def gen_walk_masks(x_0, N_particles, beta, V_xm1, V_xp1, V_x0):
    'Regner ut p minus samt p pluss og bruker uniform distribusjon til å bestemme om hver partikkel skal gå venstre, høyre eller ikke bevege seg'
    uniform_dist = np.random.uniform(0, 1, N_particles)
    go_left = uniform_dist <= (p_minus(x_0, beta, V_xm1, V_xp1, V_x0))
    go_right = uniform_dist >= 1-(p_plus(x_0, beta, V_xm1, V_xp1, V_x0))
    return go_left, go_right

@njit
def random_walk_ratchet_potential(particles, N_timesteps, T_p, N_particles, N_points, b, beta):
    #Lagrer alle partikkelstrømmer for en syklus slik at snittet kan utregnes
    particle_current_buffer = np.zeros(2*T_p)
    
    #Velger tilfeldige partikler som plottes
    plt_particle_idxs = np.random.choice(np.arange(0, N_particles), size = min(5, N_particles), replace=False)

    particle_movements = np.zeros((N_timesteps, min(5, N_particles)))
    cycle_averaged_particle_currents = np.zeros(N_timesteps // 2*T_p)
    vline_plot_points = vline_plot_points = np.zeros((N_timesteps // T_p, 2), dtype=np.int64) # Brukes i plotting senere for å indikere potensialbytte
    potential_switch_count = 0

    for j in range(N_timesteps):
        #For å regne ut partikkelstrømmen
        n_plus = 0
        n_minus = 0

        #Itererer periodisk gjennom potensialer etter T_p tidssteg
        if (j % T_p == 0) & (j != 0):
            potential_switch_count += 1
            curr_potential = potential_switch_count % 2
            vline_plot_points[j//T_p] = np.array([j, curr_potential])


        x0 = particles
        if curr_potential == 0:
            V_xm1, V_xp1, V_x0 = V_1(x0 -1), V_1(x0 + 1), V_1(x0)
        else:
            V_xm1, V_xp1, V_x0 = V_2(x0 -1), V_2(x0 + 1), V_2(x0)
        go_left, go_right = gen_walk_masks(particles, N_particles, beta, V_xm1, V_xp1, V_x0)
        movement_order = np.random.permutation(np.arange(N_particles))
        for curr_particle_idx in movement_order:
            distance_from_particle = np.copy(particles) - particles[curr_particle_idx]
            if go_left[curr_particle_idx]:
                if np.any((distance_from_particle > -b) & (distance_from_particle < 0)):
                    continue
                else:
                    #Denne if-statementen sjekker om periodisk randbetingelse er brutt
                    if curr_particle_idx <= b:
                        boundary = N_points - 1 - (b - curr_particle_idx) # -1 pga 0-indeksering
                        if np.any(distance_from_particle >= boundary):
                            continue
                
                #Denne koden kjører kun dersom ingen if-statements ovenfor er gyldig
                particles[curr_particle_idx] -= 1
                n_minus += 1

            elif go_right[curr_particle_idx]:
                if np.any((distance_from_particle < b) & (distance_from_particle > 0)):
                    continue
                if curr_particle_idx >= N_points - 1 - b:
                    boundary = b - (N_points - 1 - curr_particle_idx)
                    if np.any(np.abs(distance_from_particle) >= N_points - 1 - boundary):
                        continue

                particles[curr_particle_idx] += 1
                n_plus += 1

            #Syklisk randbetingelse for x
            if particles[curr_particle_idx] == N_points:
                particles[curr_particle_idx] = 0
            
            elif particles[curr_particle_idx] == -1:
                particles[curr_particle_idx] = N_points - 1
        
        
        normalized_particle_current = (n_plus - n_minus) / N_particles
        particle_current_buffer[j % len(particle_current_buffer)] = normalized_particle_current
        particle_movements[j] = (particles[plt_particle_idxs])

        if (j % (2*T_p) == 0) & (j != 0):
            cycle_averaged_particle_currents[N_timesteps // (2*j)] = (np.average(particle_current_buffer))

    return potential_switch_count, cycle_averaged_particle_currents, particle_movements, vline_plot_points



class ratchet_interaction_walker():
    def __init__(self, cfg : dict):
        self.cfg = cfg
        k_b = sp.constants.Boltzmann
        self.potentials = [V_1, V_2]
        self.T = cfg['T']
        self.beta = (self.T * k_b)**-1
        self.beta_k_ratio = cfg['beta_k_ratio']
        self.k = self.beta_k_ratio * (self.T * k_b)
        self.h = cfg['h']
        self.b = cfg['b']
        self.T_p = cfg['T_p']
        self.N_cycles = cfg['N_cycles']
        self.alpha = cfg['alpha']
        self.N_x = cfg['N_x']
        self.N_points = self.N_x * cfg['N_s']
        self.N_particles = cfg['N_p']
        self.N_timesteps = self.T_p * 2 * self.N_cycles

        global k
        global alpha
        global N_x
        k = self.k
        alpha = self.alpha
        N_x = self.N_x
  

    def interaction_simulator(self):
        '''Simulerer "Random walk in a ratchet potential with interactions". Returnerer alle tidssteg samt alle x-posisjoner til tilfeldig utvalgte partikler'''

        if self.N_particles == 1:
            particles = np.array([self.N_points // 2])
        else:
            particles = np.arange(0, self.N_points, self.N_points/self.N_particles, dtype=np.float64)

        self.potential_switch_count, self.cycle_averaged_particle_currents, particle_movements, self.vline_plot_points = random_walk_ratchet_potential(
            particles,
            self.N_timesteps,
            self.T_p,
            self.N_particles,
            self.N_points,
            self.b,
            self.beta
        )

        particle_movements = np.array(particle_movements).T
        self.cycle_averaged_particle_currents = np.average(self.cycle_averaged_particle_currents)
        return (np.arange(self.N_cycles * self.T_p * 2), particle_movements)
    



def oppg4a(cfg : dict):

    cfg = cfg['oppg-4a']
    walker = ratchet_interaction_walker(cfg)
    T, x_array = walker.interaction_simulator()

    for x in x_array:
        plt.plot(T, x)
    for t, potential_idx in walker.vline_plot_points:
        if potential_idx % 2 == 0:
            color = 'green'
            label = '$V_1$' if t == walker.vline_plot_points[0][0] else ''
        else:
            color = 'black'
            label = '$V_2$' if t == walker.vline_plot_points[1][0] else ''
        plt.axvline(t, color=color, linestyle='--', label = label)
    plt.legend()
    plt.show()

def rho_iterator(cfg: dict, T_p = 'default'):
    '''Finner alle cycle-averaged particle currents for oppgitt intervall av rho. Returnerer:
    cycles og liste som inneholder tuppler av (rho-verdien, liste over alle cycle-averaged particle currents for rho-verdien)'''

    plot_values = list()

    if not isinstance(T_p, str):
        cfg['T_p'] = int(T_p)
        cfg['N_c'] = 3E4 // T_p

    #Finner tilsvarende antall partikler til partikkeltettheter oppgitt i oppgaven
    N_p_min, N_p_max = np.ceil(np.array([cfg['rho_min'],cfg['rho_max']]) * cfg['N_s'] * cfg['N_x'] / cfg['b']).astype(int)
    for N_p in np.arange(N_p_min, N_p_max + 1):
        cfg['N_p'] = N_p
        rho = N_p * cfg['b'] / (cfg['N_x'] * cfg['N_s'])
        print(f'Current iteration - rho: {rho}')
        walker = ratchet_interaction_walker(cfg)
        walker.interaction_simulator()
        plot_values.append([rho, walker.cycle_averaged_particle_currents])

    return plot_values

def oppg4b(cfg : dict):
    cfg = cfg['oppg-4b']
    plot_values = rho_iterator(cfg)
    rho, avg_current = np.array(plot_values).T
    plt.plot(rho, avg_current)

    plt.title('Syklus-snittet partikkelstrømning med varierende partikkeltetthet')
    plt.xlabel('$\\rho$')
    plt.ylabel('Normalisert partikkelstrøm')
    plt.legend()
    plt.show()


def oppg4c(cfg : dict):
    T_p_vals = np.array([10, 100, 450, 1000, 10000])
    for T_p in T_p_vals:
        cycles, plot_vals = rho_iterator(cfg, T_p)