#!/usr/bin/python3
#%%
import numpy as np
import matplotlib.pyplot as plt
# %%
rdf_force = np.loadtxt('rdf_force.dat', skiprows=4)

# %%
frdf = np.loadtxt('frdf.dat')
# %%
r = rdf_force[:,1]
g1 = rdf_force[:,2]
g2 = frdf[:,1]
n = min(len(r), len(g1), len(g2))
r = r[:n]
g1 = g1[:n]
g2 = g2[:n]

# %%
plt.plot(r, g1/g2)
plt.show()
