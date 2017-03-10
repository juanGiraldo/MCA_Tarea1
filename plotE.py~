import numpy as np
import matplotlib.pyplot as plt

e=np.genfromtxt("energia.txt",delimiter=',')
tprocs=np.genfromtxt("tiempo.txt",dtype="string",usecols=0)
t=e[:,3]

#graficas de energias
plt.plot(e[:,0],t,label=r'$k=1$')
plt.plot(e[:,1],t,label=r'$k=2$')
plt.plot(e[:,3],t,label=r'$k=3$')
plt.legend(loc='best')
plt.xlabel(r'$t$')
plt.ylabel(r'$E$')
plt.title(r'$Energia\quad vs.\quad tiempo$')
plt.savefig("energias.pdf")

#grafica de tiempos de ejecucion
nprocs=np.array((1,2,4))
plt.scatter(nprocs,tprocs)
plt.xlabel("numero de procesadores")
plt.ylabel("tiempo de ejecucion")
plt.title("#procesadores vs. tiempo")
plt.savefig("tiempos.pdf")
