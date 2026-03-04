import numpy as np
from numba import njit

@njit
def p_plus(beta, V_xm1, V_xp1, V_x0):
    denonimator = 1 + np.exp(-beta*(V_xm1 - V_xp1)) + np.exp(-beta*(V_x0 - V_xp1))
    return 1/denonimator

@njit
def p_zero(beta, V_xm1, V_xp1, V_x0):
    denonimator = 1 + np.exp(-beta*(V_xm1 - V_x0)) + np.exp(-beta*(V_xp1 - V_x0))
    return 1/denonimator

@njit
def p_minus(beta, V_xm1, V_xp1, V_x0):
    denonimator = 1 + np.exp(-beta*(V_xp1 - V_xm1)) + np.exp(-beta*(V_x0 - V_xm1))
    return 1/denonimator