import scipy as sp
import numpy as np

from utilities.probs import p_minus, p_plus



alpha = 0.8
T = 273.15 + 37
beta = sp.constants.Boltzmann *T
N_x = 100
N_points = 2*N_x
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
        return k * x * (alpha*N_x)
    if(-(1-alpha)*N_x < x <= 0):
        return -k*x*((1-alpha)*N_x)

def V_2(x):
    return k

def random_walk_particle(x : np.array, walker : int, V):
    uniform = np.random.uniform()
    step_dir = 0
    if(uniform > p_plus(x[walker], beta, k, V)):
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


def random_walk_cycle(positions,occupied):
    avg_stream = 0
    for timestep in range(2*T_p):
        potential = V_2 #antar V1 med mindre vi er i en oddetalls del av syklusen
        if (timestep >= T_p):
            potential = V_1
        occupied, stream = random_walk_system(positions, occupied, potential)
        avg_stream += stream

    avg_stream /= 2*T_p

    return occupied, avg_stream

def oppgave_3_a(num_cycles = 10):
    positions = np.linspace(0, h*(N_points-1), N_points)
    occupied = np.ones((N_points),dtype=np.int16)*N_particles/N_points
    streams = np.zeros(num_cycles)
    for i in range(num_cycles):
        occupied, streams[i] = random_walk_cycle(positions, occupied)
        print(f"cycle: {i+1}, average stream {streams[i]}")

    