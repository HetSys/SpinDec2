#!/usr/bin/env python3

import netCDF4 as NC
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from scipy.ndimage import uniform_filter1d
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse

plt.rcParams["figure.figsize"] = [8, 8]
plt.rcParams["figure.autolayout"] = True
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14



data = NC.Dataset('CH_output.nc', "r", format="NETCDF4")

c = np.array([data.variables['c'][:]][0])[0]
kappa = data.kappa
Nx = data.Nx

dx = 1/Nx

x = np.linspace(0, 1, Nx)

print(x)

# fig, ax = plt.subplots()
# pos1 = ax.imshow(c, vmin=0, vmax=1, origin = 'lower')
# div = make_axes_locatable(ax)
# cax = div.append_axes('right', '5%', '5%')
# fig.colorbar(pos1, cax = cax)
      
# ax.set_xlabel(r'x')
# ax.set_ylabel(r"$y$")
# ax.set_title(r'$c(x,y)$')


# plt.show()


def profile(x, x0, kappa, cA = 1, cB = 0):

    phi = (cA+cB)/2 + (cA-cB)/2 * np.tanh((x-x0) / np.sqrt(2*kappa))

    return phi


pos = np.where(c[64,:] > 1e-3)
x0 = x[pos][6]

dc = np.diff(c[64,:])
print(dc)

plt.plot(dc)
plt.show()

print(dc)

print(x0)


# fig, ax = plt.subplots()

# ax.scatter(x, c[64, :], color = 'C1')
# ax.plot(x, profile(x, x0, kappa=kappa), color = 'C2')
# ax.legend()
# plt.xlim(0.2, 0.6)

# plt.show()