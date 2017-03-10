import numpy as np
import atplotlib.pyplot as plt

et=np.genfromtxt("energias.txt",delimiter=',')
t=e[:,3]
plt.plot(e[:,0],t,label=r'$k=1$')
plt.plot(e[:,1],t,label=r'$k=2$')
plt.plot(e[:,3],t,label=r'$k=3$')
plt.legend(loc='best')
plt.xaxis(r'$t$')
plt.yaxis(r'$E$')
plt.title(r'$Energia\quad vs.\quad tiempo$')
