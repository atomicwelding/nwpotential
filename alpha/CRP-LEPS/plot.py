import matplotlib.pyplot as plt
import numpy as np

step = 35

high_sym_sites = ["top", "bridge", "hollow"]
low_sym_sites =  ["topbridge", "tophollow", "bridgehollow"]

colors = ["#f2bc46", "#bc3754", "#460aad"]
markers = ["^", "s", "p"]

# iterate over each sites
plt.figure(figsize=(8,8))
for idx, site in enumerate(low_sym_sites):

    # loading data
    file_crp = site + "crp.txt"
    file_leps = site + "leps.txt"

    crp = np.loadtxt(file_crp)
    leps = np.loadtxt(file_leps, skiprows=1)

    # computing error
    errors = leps[:,1] - crp[:,1]
    mean_error = np.mean(errors)
    absolute_deviation = np.mean(np.abs(errors))

    print(f"{site}")
    print(f"mean error (LEPS - CRP): {mean_error:.6f} eV")
    print(f"absolute deviation : {absolute_deviation:.6f} eV")

    # plotting
    plt.plot(leps[:,0], leps[:,1], color=colors[idx])
    plt.plot(crp[:,0][::step], crp[:,1][::step],
              linestyle='dotted',
              color=colors[idx],
              marker=markers[idx],
              markersize=8,
              label=site)

# fake lines to add to the legend
plt.plot([-1,-0.9], [3,4], 'k*', linestyle="dotted",label="CRP")
plt.plot([-1,-0.9], [3,4], 'k', label="LEPS")

# ending plot configuration and show
plt.xticks(fontsize=16) 
plt.yticks(fontsize=16)
plt.legend(fontsize=18)
plt.xlabel(r"$z_a$ ($\AA$)", fontsize=18)
plt.ylabel(r"$V(z_a)$ (eV)", fontsize=18)
plt.ylim(top=2, bottom=-8)
plt.xlim(left=0., right=6.)
plt.grid()
plt.show()
