import numpy as np
import matplotlib.pyplot as plt
file=np.genfromtxt("data.txt")
plt.hist(file,bins='auto')
plt.title("decaimiento exponencial")
plt.savefig("histograma_decaimiento.pdf")
