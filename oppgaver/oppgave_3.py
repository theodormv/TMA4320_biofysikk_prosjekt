import scipy as sp
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from utilities.probs import p_minus, p_plus


import matplotlib.animation as animation
import functools


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



def V_1(x):
    #normaliserer x til periodisiteten
    while x > alpha*N_x:
        x = x - alpha*N_x
    while x < -(1-alpha)*N_x:
        x = x + alpha*N_x

    #delt funskjnons utrykk
    if(0 < x <= (alpha*N_x)):
        return k * x / (alpha*N_x)
    if(-(1-alpha)*N_x < x <= 0):
        return -k*x/((1-alpha)*N_x)
    
def V_1_vectorized(x):
    #normaliserer x til periodisiteten
    V_1 = np.zeros(np.shape(x))

    L = -(1-alpha)*N_x
    positions = L + (x-L)%N_x
    
    
    V_1 += (positions > 0) * positions*k / (alpha*N_x)
    V_1 += -1*(positions <= 0)*positions*k/ ((1-alpha)*N_x)
    
    return V_1
    

def V_2(x):
    return k



def random_walk_particle(x : np.array, walker : int, V):
    uniform = np.random.uniform()
    step_dir = 0
    if(uniform > 1 - p_plus(x[walker], beta, k, V)):
        step_dir = 1
    if(uniform  < p_minus(x[walker], beta, k, V)):
        step_dir = -1

    new_ind = walker
    if (step_dir == 1 and walker == (N_points-1)):
        new_ind = 0
    else:
        new_ind = walker + step_dir
    
    return new_ind, step_dir


def random_walk_system(positions, occupied, potential):
    n_plus = 0
    n_minus = 0

    new_occupied = np.zeros(np.shape(occupied), dtype=np.int16)
    for index, particles  in enumerate(occupied):
        for _ in range(int(particles)):
            new_index, step_direction = random_walk_particle(positions ,index, potential)
            new_occupied[new_index] += 1
            if(step_direction == 1):
                n_plus += 1
                continue
            if(step_direction == -1):
                n_minus += 1
                continue
    
    stream = (n_plus - n_minus)/N_particles
    
    return new_occupied, stream

def random_walk_system_vectorized(particles, positions, potential):
    p = np.copy(particles)
    movements = np.random.uniform(size=N_particles)
    particle_probs = np.array((p_minus(positions[particles], beta, potential), p_plus(positions[particles], beta, potential)))

    move_left = np.less_equal(movements, particle_probs[0])*1
    move_right = np.greater_equal(movements, 1 - particle_probs[1])*1
    n_plus = np.sum(move_right)
    n_minus = np.sum(move_left)

    particles += move_right - move_left
    
    #if(np.linalg.norm(np.equal(particles, N_points))):

        #print(np.linalg.norm(np.equal(particles, N_points)))

    particles -= np.equal(particles, N_points)*N_points
    particles += np.equal(particles, -1)*N_points


    return particles, (n_plus-n_minus) / N_particles
    


def random_walk_cycle(positions,occupied):
    avg_stream = 0
    for timestep in range(2*T_p):
        if (timestep == 0 or timestep == T_p):
            plt.hist(positions, occupied)

        potential = V_2 #antar V1 med mindre vi er i en oddetalls del av syklusen
        if (timestep >= T_p):
            potential = V_1
        occupied, stream = random_walk_system(positions, occupied, potential)
        avg_stream += stream

    avg_stream /= 2*T_p

    return occupied, avg_stream

def random_walk_cycle_vectorized(positions,particles, particle_archive):
    avg_stream = 0
    
    for timestep in range(2*T_p):
        
        #antar V1 med mindre vi er i en oddetalls del av syklusen
        V = V_2
        if timestep > T_p:
            V = V_1_vectorized

        particles, stream = random_walk_system_vectorized(particles, positions, V)

        if (timestep % int(T_p/25) == 0):
            particle_archive.append(np.copy(particles))

        avg_stream += stream

    avg_stream /= 2*T_p

    return particles, avg_stream

def oppgave_3_a(num_cycles = 10):
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

    #set up kode for animasjon
    HIST_BINS = np.linspace(0, N_points*h - 1, N_points)

    particle_archive = [particles]
    

    # hoved syklus
    for cycle in range(num_cycles):
        particles, streams[cycle] = random_walk_cycle_vectorized(positions, particles, particle_archive)
        print(f"cycle: {cycle}, average stream {streams[cycle]}")
    
    """Primary animation code"""

    print(len(particle_archive))

    def animate(frame_number, bar_container):
        # Simulate new data coming in.
        n, _ = np.histogram(particle_archive[frame_number], HIST_BINS)
        for count, rect in zip(n, bar_container.patches):
            rect.set_height(count)

        return bar_container.patches
    

    fig, ax = plt.subplots()
    _, _, bar_container = ax.hist(particle_archive[0], HIST_BINS, lw=1,
                                ec="yellow", fc="green", alpha=0.5)


    anim = functools.partial(animate, bar_container=bar_container)
    ani = animation.FuncAnimation(fig, anim, 400, repeat=False, blit=True)
    plt.show()


def oppgave_3_b(num_values = 50):
    global N_particles
    N_particles = 4*N_x

    global alpha
    alpha = 0.8

    particles_start = np.zeros(N_particles, dtype=np.int16)
    particles_start[::2] = h*N_x
    positions = np.linspace(0, h*N_points, N_points)
    streams = np.zeros(num_values)
    i = 0

    global T_p
    Tp_range = np.linspace(1, 1001, num_values, dtype=np.int16)
    for T in tqdm(Tp_range):
        T_p = T
        particles = np.copy(particles_start)
        _, streams[i] = random_walk_cycle_vectorized(positions, particles)
        i+=1



    plt.plot(Tp_range, streams)
    #plt.plot(positions, V_1_vectorized(positions))
    plt.grid()
    plt.show()



def analytical_avg_current_alpha(a):

    multiplier = N_x/2 * np.sqrt(3/T_p)
    erfc_left = sp.special.erfc(a*multiplier)
    erfc_right = sp.special.erfc((1-a)*multiplier)

    return N_x/(4*T_p) * (erfc_left - erfc_right)


def analytical_avg_current_T_p(Tp):

    multiplier = N_x/2 * np.sqrt(3/Tp)
    erfc_left = sp.special.erfc(alpha*multiplier)
    erfc_right = sp.special.erfc((1-alpha)*multiplier)

    return N_x/(4*T_p) * (erfc_left - erfc_right)


def oppgave_3_c(num_values = 50):
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
        _, streams[i] = random_walk_cycle_vectorized(positions, particles)
        i+=1

    plt.plot(alpha_range, analytical_avg_current_alpha(alpha_range), label=r"Analytical $J_\text{avg}$")
    plt.plot(alpha_range, streams, label=r"Simulated $J_\text{avg}$")
    plt.legend()
    plt.grid()
    plt.show()
    

def oppgave_3_d(num_values = 50):
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
            _, streams[i] = random_walk_cycle_vectorized(positions, particles)

        
        plt.plot(alpha_range, streams, label=r"$\beta k = $" + str(kbeta[j]))
    plt.plot(alpha_range, analytical_avg_current_alpha(alpha_range), label=r"Analytical $J_\text{avg}$")
    plt.legend()
    plt.grid()
    plt.show()
    

def oppgave_3_e(num_values = 20):

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
        _, streams[i] = random_walk_cycle_vectorized(positions, particles)



    plt.plot(Tp_range, streams, label=r"simulatet $J_\text{avg}$") 
    plt.plot(Tp_range, analytical_avg_current_T_p(Tp_range), label=r"analytical $J_\text{avg}$")
    plt.legend()
    plt.grid()
    plt.show()