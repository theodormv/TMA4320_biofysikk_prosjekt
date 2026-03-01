from oppgaver.oppgave_3 import oppgave_3_a, oppgave_3_b
from oppgaver.oppgave_4_fast import oppg4a, oppg4b, oppg4c, V_1, ratchet_interaction_walker
from oppgaver.oppgave_2 import oppg2
import matplotlib.pyplot as plt
import numpy as np

cfg = {"oppg-4a":{
  "alpha": 0.2,
  "T": 310.15,
  "N_x": 20,
  "T_p": 40,
  "N_cycles": 5,
  "N_s": 4,
  "h": 1,
  "beta_k_ratio": 1000,
  "N_p": 5,
  "b": 2}, # Partikkelstørrelse

"oppg-4b":{
  "alpha": 0.2,
  "T": 310.15,
  "N_x": 100,
  "T_p": 300,
  "N_cycles": 100,
  "N_s": 10,
  "h": 1,
  "beta_k_ratio": 1000,
  "N_p": 34,
  "b": 20, # Partikkelstørrelse
  "rho_min": 0.01,
  "rho_max": 1}}


def main():
   oppg4c(cfg)
   #oppgave_3_a()

if __name__ == "__main__":
   main()