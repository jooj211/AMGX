import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-1,1,9)
xx = np.linspace(-1,1,9)
y = 1 + xx*xx

plt.plot(xx,y,label='exata')

y_a = np.loadtxt('../build/plotSol.csv')
plt.plot(x,y_a,label='aproximada')


plt.legend()
plt.show()
#plt.savefig('../build/normaL2.png')