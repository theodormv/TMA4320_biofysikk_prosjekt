import numpy as np

def p_plus(x_0, beta, k, V):
    denonimator = 1 + np.exp(-beta*(V(x_0 - 1, k) - V(x_0 + 1, k))) + np.exp(-beta*(V(x_0, k) - V(x_0 + 1, k)))
    return 1/denonimator

def p_zero(x_0, beta, k, V):
    denonimator = 1 + np.exp(-beta*(V(x_0 - 1, k) - V(x_0, k))) + np.exp(-beta*(V(x_0 + 1, k) - V(x_0, k)))
    return 1/denonimator

def p_minus(x_0, beta, k, V):
    denonimator = 1 + np.exp(-beta*(V(x_0 + 1, k) - V(x_0 - 1, k))) + np.exp(-beta*(V(x_0, k) - V(x_0 - 1, k)))
    return 1/denonimator