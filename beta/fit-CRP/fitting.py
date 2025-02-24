import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit

# utils
Ntot = 1000
zmin = 0.
zmax = 6.

def morse_potential(z, D, a0_rep, a1_rep, a0_att, a1_att, req, a, b):
    f = 0.5 * ( 1 + np.tanh(a*z + b))
    alpha = (1-f)*(a0_rep + a1_rep*z) + f*(a0_att + a1_att*z)
    return D * ( np.exp(-2 * alpha*(z - req)) - 2 * np.exp(- alpha * (z - req)) )

def printp(name, p, var):
    print(f"{name}={p:.6f}")
    print(f"{name}var={var:.6f}")

# fitting
FILE = "tophollowcrp.txt"
data = np.loadtxt(FILE)


## 220 for z = 1.32
zs = data[:,0]
V = data[:,1]

# boundaries
lower_bounds = [-np.inf, -np.inf, -np.inf, -np.inf, -np.inf, 0, 0, -5]
upper_bounds = [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,  5, 0]

## proper fitting
popt, pcov = curve_fit(morse_potential, zs, V,
                       bounds=(lower_bounds, upper_bounds))
D, a0_rep, a1_rep, a0_att, a1_att, req, a, b = popt

print(f"D       = {D:.6f}")
print(f"a0_rep  = {a0_rep:.6f}")
print(f"a1_rep  = {a1_rep:.6f}")
print(f"a0_att  = {a0_att:.6f}")
print(f"a1_att  = {a1_att:.6f}")
print(f"req     = {req:.6f}")
print(f"a       = {a:.6f}")
print(f"b       = {b:.6f}")

## plotting
plt.figure(figsize=(8,8))
plt.plot(zs, V, 'k', label="CRP")
plt.plot(zs, morse_potential(zs, D, a0_rep, a1_rep, a0_att, a1_att, req, a, b),
         'r', label="morse fit")
plt.legend(fontsize=14)
plt.xlabel(r"$z_a$  ($\AA$)", fontsize=16)
plt.ylabel(r"$V(z_a)$ (eV)", fontsize=16)
plt.title(FILE)
plt.grid()
plt.xlim(left=zs[0], right=zs[-1])
plt.show()
