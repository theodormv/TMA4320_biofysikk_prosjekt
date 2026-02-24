from oppgaver.oppgave_3 import oppgave_3_a
from oppgaver.oppgave_4 import ratchet_interaction_walker, cfg
import matplotlib.pyplot as plt

def main():



   #Oppg 4
   walker = ratchet_interaction_walker(cfg)
   T, x_array = walker.interaction_simulator()

   for x in x_array:
      plt.plot(T, x)
   for i in range(len(x)):
      if i % 40 == 0:
         if i % 80 == 0:
            color = 'green'
         else:
            color = 'black'
         plt.axvline(i, color=color, linestyle='--')
   plt.show()


if __name__ == "__main__":
    main()