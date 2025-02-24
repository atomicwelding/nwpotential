import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  

filename = '2dcut.txt'

data = np.loadtxt(filename, skiprows=1)

x = data[:, 0]
y = data[:, 1]
V = data[:, 2]

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_trisurf(x, y, V, cmap='magma', linewidth=0.2)

ax.set_xlabel(r'$x$ (in units of $\delta$)', fontsize=18)
ax.set_ylabel(r'$y$ (in units of $\delta$)', fontsize=18)
ax.set_zlabel(r'Potential energy in eV', fontsize=18)

plt.xticks(fontsize=12) 
plt.yticks(fontsize=12)
ax.tick_params(axis='z', labelsize=12)

plt.tight_layout()
plt.show()
