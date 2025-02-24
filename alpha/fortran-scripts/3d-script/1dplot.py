import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("1dcut.txt", skiprows=1)
zs = data[:,0]
V  = data[:,1]

plt.figure(figsize=(8,8))
plt.plot(zs, V, 'k')
plt.xlabel(r'$z$ (in angstroms)', fontsize=18)
plt.ylabel('Potential energy in eV', fontsize=18)
plt.xticks(fontsize=14)
plt.xlim(left=zs[0], right=zs[-1])
plt.yticks(fontsize=14)
plt.grid()
plt.show()

