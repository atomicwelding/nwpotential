import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit

# utils
Ntot = 10_000
zmin = 0.
zmax = 6.

def morse_potential(z, D = 1, alpha = 1, req = 1):
    return D * ( np.exp(-2 * alpha*(z - req)) - 2 * np.exp(- alpha * (z - req)) )

def printp(name, p, var):
    print(f"{name}={p:.6f}")
    print(f"{name}var={var:.6f}")

# fitting
FILE = "top.dat"
data = np.loadtxt(FILE)

zs = data[2100:,0]
V = data[2100:,1]

## guessing parameters
min_idx = np.argmin(V)
Dguess = V[min_idx]
reqguess = zs[min_idx]
alphaguess = 1/reqguess

## proper fitting
popt, pcov = curve_fit(morse_potential, zs, V,
                       p0=[Dguess, alphaguess, reqguess],
                       method='lm')
D, alpha, req = popt
Dvar, alphavar, reqvar = pcov[0,0], pcov[1,1], pcov[2,2]

printp("D", D, Dvar)
printp("alpha", alpha, alphavar)
printp("req", req, reqvar)

## plotting
plt.figure(figsize=(8,8))
plt.plot(zs, V, 'k', label="CRP")
plt.plot(zs, morse_potential(zs, D, alpha, req), 'r', label="morse fit")
plt.legend(fontsize=14)
plt.xlabel(r"$z_a$  ($\AA$)", fontsize=16)
plt.ylabel(r"$V(z_a)$ (eV)", fontsize=16)
plt.title(FILE)
plt.grid()
plt.xlim(left=zs[0], right=zs[-1])
plt.show()
