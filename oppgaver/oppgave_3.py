import scipy as sp
import numpy as np



alpha = 0.8
T = 273.15 + 37
beta = sp.constants.Boltzmann *T
N_x = 100
N_points = 2*N_x
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

def random_walk_index_in_rachet(x : np.array, occupied : np.array, walker : int, V : function)
    uniform = np.random.uniform()
    step_dir = 0
    if(uniform > p_plus(x[walker], beta, k, V(x[walker]))):
        step_dir = 1
    if(uniform  < p_minus(x[walker], beta, k, V(x[walker]))):
        step_dir = -1

    new_ind = walker
    if (step_dir == 1 and walker == (N_points-1)):
        new_ind = 0
    else:
        new_ind = walker + step_dir
    
    return new_ind


    