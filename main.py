from oppgaver.oppgave_3 import oppgave_3_a
from oppgaver.oppgave_4 import ratchet_interaction_walker, cfg
import matplotlib.pyplot as plt

def main():



   #Oppg 4
   walker = ratchet_interaction_walker(cfg)
   x, y = walker.interaction_simulator()

   plt.plot(x, y)
   for i in x:
       if i % 40 == 0:
           plt.axvline(i, color='black', linestyle='--')
   plt.show()


if __name__ == "__main__":
    main()