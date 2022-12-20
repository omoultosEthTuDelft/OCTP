#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# find all files in the current directory that start with "frdf" and end in "".dat"

frdf_files = sorted(Path(".").glob("frdf*.dat"))
rdf_force_files = sorted(Path(".").glob("rdf_force*.dat"))
lmp_force_files = sorted(Path(".").glob("rdf_lmp*.dat"))

for file in frdf_files:
    print(file)
    data = np.loadtxt(file)
    plt.plot(data[:,0], data[:,1], label=file.stem)

for file in rdf_force_files:
    print(file)
    data = np.loadtxt(file, skiprows=4)
    # plt.plot(data[:,1], data[:,2]*(10**3-1)/10**3, label=file.stem+"_corrected")
    plt.plot(data[:,1], data[:,2], label=file.stem+"lmp")
    plt.plot(data[:,1], data[:,4], label=file.stem+"force")

for file in lmp_force_files:
    print(file)
    data = np.loadtxt(file, skiprows=4)
    # plt.plot(data[:,1], data[:,2], label=file.stem)
    plt.plot(data[:,1], data[:,2]*(12**3-1)/12**3, label=file.stem)

# %%
try:
    rdf_ext = np.loadtxt('rdf_ext.dat')
    with open('rdf_ext.dat') as f:
        # This is to read the column names
        header = f.readline().strip().split()[1:]

    # plt.figure()
    for i in range(1, rdf_ext.shape[-1],2):
        plt.plot(rdf_ext[:,0], rdf_ext[:,i], label=header[i])
        plt.xlabel(header[0])

    # %%
    for i in range(2, rdf_ext.shape[-1],2):
        plt.plot(rdf_ext[:,0], rdf_ext[:,i], label=header[i])
        plt.xlabel(header[0])
except FileNotFoundError:
    pass

plt.grid()
# add small grid with thinner lines
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.minorticks_on()

plt.legend()
plt.show()    