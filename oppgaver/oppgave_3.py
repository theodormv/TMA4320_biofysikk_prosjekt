import scipy as sp
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from utilities.probs import p_minus, p_plus


alpha = 0.8
T = 273.15 + 37
beta = sp.constants.Boltzmann *T
N_x = 100
periods = 2
N_points = periods*N_x
N_particles = 12*N_x
k = 1000/beta

T_p = 500
h = 1

def V_1_vectorized(x, cfg):
    #normaliserer x til periodisiteten
    V_1 = np.zeros(np.shape(x))

    L = -(1-alpha)*N_x
    positions = L + (x-L)%N_x

    
    V_1 += (positions > 0) * positions*k / (alpha*N_x)
    V_1 += -1*(positions <= 0)*positions*k/ ((1-alpha)*N_x)
    
    return V_1
    

def V_2(x):
    return k



def random_walk_system_vectorized(particles, positions, potential, cfg):

    movements = np.random.uniform(size=N_particles)
    particle_probs = np.array((p_minus(positions[particles], beta, potential), p_plus(positions[particles], beta, potential)))

    move_left = np.less_equal(movements, particle_probs[0])*1
    move_right = np.greater_equal(movements, 1 - particle_probs[1])*1
    n_plus = np.sum(move_right)
    n_minus = np.sum(move_left)

    particles += move_right - move_left

    
    particles -= np.equal(particles, N_points)*N_points
    particles += np.equal(particles, -1)*N_points

    for i in range(len(particles)):
        if(particles[i] == 200 or particles[i]==-1):
            print(f"index: {i} == {particles[i]}")

    return particles, (n_plus-n_minus) / N_particles
    


def random_walk_cycle_vectorized(positions,particles, cfg):
    avg_stream = 0
    
    for timestep in range(2*T_p):
        
        #antar V1 med mindre vi er i en oddetalls del av syklusen
        V = V_2
        if timestep > T_p:
            V = V_1_vectorized

        

        particles, stream = random_walk_system_vectorized(particles, positions, V, cfg)
        avg_stream += stream

    avg_stream /= 2*T_p

    return particles, avg_stream


def analytical_avg_current(a, Tp, cfg):

    multiplier = N_x/2 * np.sqrt(3/Tp)
    erfc_left = sp.special.erfc(a*multiplier)
    erfc_right = sp.special.erfc((1-a)*multiplier)

    return N_x/(4*Tp) * (erfc_left - erfc_right)


def oppgave_3_a(cfg , num_cycles = 10):
    global N_x
    N_x = 100
    
    global N_particles
    N_particles = 12*N_x

    global alpha
    alpha = 0.8

    global T_p
    T_p = 500

    #initialliserer partikellene jevnt utover 
    particles = np.array((np.arange(0, N_points, 1, dtype=np.int16),))
    one = np.ones((int(N_particles/N_points),1), dtype = np.int16)
    particles = np.ndarray.flatten(one @ particles)

    positions = np.linspace(0, h*N_points-1, N_points)
    streams = np.zeros(num_cycles)

    for cycle in range(num_cycles):
        particles, streams[cycle] = random_walk_cycle_vectorized(positions, particles)
        #plt.hist(particles, bins = 20)
        print(f"cycle: {cycle}, average stream {streams[cycle]}")
        #plt.show()
   

def oppgave_3_b(cfg , num_values = 50):
    global N_particles
    N_particles = 4*N_x

    global alpha
    alpha = 0.8

    particles_start = np.zeros(N_particles, dtype=np.int16)
    particles_start[::2] = h*N_x
    positions = np.linspace(0, h*N_points, N_points)
    streams = np.zeros(num_values)
    

    global T_p
    Tp_range = np.linspace(1, 1001, num_values, dtype=np.int16)
    for T in tqdm(Tp_range):
        T_p = T
        particles = np.copy(particles_start)
        _, streams[i] = random_walk_cycle_vectorized(positions, particles, cfg)
        i+=1



    plt.plot(Tp_range, streams)
    #plt.plot(positions, V_1_vectorized(positions))
    plt.grid()
    plt.title(f"gjennomsnittslig simulert strømnig mned variert $T_p$")
    plt.xlabel(f"$T_p$")
    plt.ylabel(r"$J_\text{avg}$")
    plt.show()





def oppgave_3_c(cfg , num_values = 50):
    global T_p
    T_p = 500
    global N_particles
    N_particles = 12*N_x
    
    particles_start = np.zeros(N_particles, dtype=np.int16)
    particles_start[::2] = h*N_x
    positions = np.linspace(0, h*N_points, N_points)
    streams = np.zeros(num_values)
    i = 0

    global alpha
    alpha_range = np.linspace(1e-5,1-1e-5,num_values) #Velger små verdier som approximerer de egentlige grensene for å unngå overflow 

    for a in tqdm(alpha_range):
        alpha = a
        particles = np.copy(particles_start)
        _, streams[i] = random_walk_cycle_vectorized(positions, particles, cfg)
        i+=1

    plt.plot(alpha_range, analytical_avg_current(alpha_range, cfg.T_p, cfg), label=r"Analytical $J_\text{avg}$")
    plt.plot(alpha_range, streams, label=r"Simulated $J_\text{avg}$")
    plt.legend()
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"$J_\text{avg}$")
    plt.title("analytisk og simulert snitt strømning")
    plt.grid()
    plt.show()
    

def oppgave_3_d(cfg ,num_values = 50):
    global T_p
    T_p = 500
    global N_particles
    N_particles = 12*N_x
    
    particles_start = np.zeros(N_particles, dtype=np.int16)
    particles_start[::2] = h*N_x
    positions = np.linspace(0, h*N_points, N_points)
    streams = np.zeros(num_values)
    

    global alpha
    alpha_range = np.linspace(1e-5,1-1e-5,num_values) #Velger små verdier som approximerer de egentlige grensene for å unngå overflow 

    global k
    kbeta = np.array((0.01,1,2,3,5,10))
    k_range = kbeta/beta
    np.set_printoptions(precision=5)
    for j,e in enumerate(tqdm(k_range)):
        k = e
        for i,a in enumerate(alpha_range):
            alpha = a
            particles = np.copy(particles_start)
            _, streams[i] = random_walk_cycle_vectorized(positions, particles, cfg)

        
        plt.plot(alpha_range, streams, label=r"$\beta k = $" + str(kbeta[j]))
    plt.plot(alpha_range, analytical_avg_current_alpha(alpha_range), label=r"Analytical $J_\text{avg}$")
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"$J_\text{avg}$")
    plt.title("sammenlikning av analytisk og simulert strømning for forskjellige forhold ")
    plt.legend()
    plt.grid()
    plt.show()
    

def oppgave_3_e(cfg ,num_values = 20):

    global N_x
    N_x = 100

    global N_points
    N_points = periods*N_x

    global N_particles
    N_particles = 40*N_x

    global alpha
    alpha = 0.8


    particles_start = np.zeros(N_particles, dtype=np.int16)
    particles_start[::2] = h*N_x
    positions = np.linspace(0, h*N_points, N_points)
    streams = np.zeros(num_values)

    global T_p
    Tp_range = np.linspace(80, 1500, num_values, dtype=np.int16)
    for i, T in enumerate(tqdm(Tp_range)):
        T_p = T
        particles = np.copy(particles_start)
        _, streams[i] = random_walk_cycle_vectorized(positions, particles, cfg)



    plt.plot(Tp_range, streams, label=r"simulatet $J_\text{avg}$") 
    plt.plot(Tp_range, analytical_avg_current(cfg.alpha, Tp_range, cfg), label=r"analytical $J_\text{avg}$")
    plt.legend()
    plt.grid()

    plt.show()