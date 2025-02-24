import numpy as np
import matplotlib.pyplot as plt


def morse_potential(z, D = 1, alpha = 1, req = 1):
    return D * ( np.exp(-2 * alpha*(z - req)) - 2 * np.exp(- alpha * (z - req)) )


z = np.linspace(0., 6.,  num=200)
V = morse_potential(z)


plt.figure(figsize=(8,8))
plt.plot(z,V,'k')

plt.xlabel(r"$Z$  ($\AA$)", fontsize=16)
plt.ylabel(r"$V(Z)$ (eV)", fontsize=16)

plt.grid()
plt.xlim(left=z[0], right=z[-1])
plt.show()
