import numpy as np
import matplotlib.pyplot as plt


FILE = "2dcut.txt"
data = np.loadtxt(FILE, skiprows = 1)


x = data[:,0]
y = data[:,1]
V = data[:,2]

fig, ax = plt.subplots(figsize=(8,8))
contour = ax.tricontour(x,y,V, cmap='magma')
ax.clabel(contour, inline=True, fontsize=15)
plt.xlabel(r"$x$ (in units of $\delta$)", fontsize=18)
plt.ylabel(r"$y$ (in units of $\delta$)", fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.grid()
plt.show()
