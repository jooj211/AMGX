import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

size = int(sys.argv[1])
figure_name = sys.argv[2] + "_" + sys.argv[3] + ".png"

x = np.linspace(-1, 1, size)
xx = np.linspace(-1, 1, size)
y = 1 + xx * xx

plt.plot(xx, y, label='exata')

y_a = np.loadtxt('../build/plotSol.csv')
plt.plot(x, y_a, label='aproximada')

plt.legend()
plt.savefig("fig/" + figure_name)
