#!/usr/bin/env python3

import netCDF4 as NC
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from scipy.ndimage import uniform_filter1d
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse


from scipy.optimize import curve_fit
plt.rcParams["figure.figsize"] = [8, 8]
plt.rcParams["figure.autolayout"] = True
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 16



data = NC.Dataset('CH_output.nc', "r", format="NETCDF4")

c = np.array([data.variables['c'][:]][0])[0]
kappa = data.kappa
Nx = data.Nx

dx = 1/Nx

x = np.linspace(0, 1, Nx)



fig, ax = plt.subplots()
pos1 = ax.imshow(c, vmin=0, vmax=1, origin = 'lower')
div = make_axes_locatable(ax)
cax = div.append_axes('right', '5%', '5%')
fig.colorbar(pos1, cax = cax)
      
ax.set_xlabel(r'x')
ax.set_ylabel(r"$y$")
ax.set_title(r'$c(x,y)$')



plt.show()


def profile(x, x0, kappa, cA = 1, cB = 0):

    phi = (cA+cB)/2 + (cA-cB)/2 * np.tanh((x-x0) / np.sqrt(2*kappa))

    return phi


# pos = np.where(c[250,:] < 0.55)
# x0 = x[pos][-1]

# # print(pos)


# fig, ax = plt.subplots()

# ax.scatter(x, c[250, :], color = 'C0', label = r'SpinDec2 Result')
# ax.plot(x, profile(x, x0 = x0, kappa=kappa), color = 'C1', label = r'$\frac{c_A+c_B}{2} + \frac{c_A-c_B}{2} \,\, \tanh\left(\frac{(x-x_0)}{\sqrt{2 \kappa}}\right)$')
# ax.legend(loc='upper left', frameon = False)
# ax.text(x = 0.56, y = 0.85, s = r'$c_A = 1.0, c_B = 1.0$', fontsize = 14)
# ax.text(x = 0.56, y = 0.8, s = r'$\kappa = $' + str(kappa), fontsize = 14)
# ax.text(x = 0.56, y = 0.75, s = r'$x_0 = $' + str(np.round(x0, 4)), fontsize = 14)
# # ax.axvspan(xmin=0.66, xmax = 0.7, alpha = 0.4)
# plt.xlim(0.55, 0.8)

# plt.show()




# fig, ax = plt.subplots(1, 3, figsize=(20,20))


# locs = [0, Nx//2, 499]
# col = ['C0', 'C1', 'C2']

# for j in range(3):
#     pos = np.where(c[locs[j],:] < 0.55)
#     x0 = x[pos][-1]
#     ax[j].scatter(x, c[locs[j], :], color = col[j],  label = r'SpinDec2 Result')
#     ax[j].plot(x, profile(x, x0 = x0, kappa=kappa), linestyle = '--', color = 'black', label = r'$A + B \,\, \tanh\left(\frac{(x-x_0)}{\sqrt{2 \kappa}}\right)$')
#     ax[j].set_xlim(x0*0.75, 1.25*x0)
#     ax[j].legend(loc = 'upper left', frameon = False)
#     ax[j].text(x = x0*0.76, y = 0.85, s = r'$A = \frac{c_A + c_B}{2} \,\, , B = \frac{c_A-c_B}{2}$', fontsize = 15)
#     ax[j].text(x = x0*0.76, y = 0.80, s = r'$c_A = 1.0, c_B = 1.0$', fontsize = 13)
#     ax[j].text(x = x0*0.76, y = 0.75, s = r'$\kappa = $' + str(kappa), fontsize = 13)
#     ax[j].text(x = x0*0.76, y = 0.7, s = r'$x_0 = $' + str(np.round(x0, 4)), fontsize = 13)
#     ax[j].text(x = x0*0.76, y = 0.65, s = r'$N_y = $' + str(np.round(locs[j], 4)), fontsize = 13)
#     ax[j].set_xlabel(r'$x$')
#     ax[j].set_ylabel(r'$c(x)$')
# plt.show()


np.save('kappa_'+str(kappa)+'.npy', c[Nx//2, :])

kappas = [1e-05, 3e-05, 5e-05, 8e-05, 0.0001, 0.0003, 0.0005, 0.0008, 0.001, 0.003, 0.005, 0.008, 0.01]

sig_data = np.zeros((len(kappas),330))

x = np.linspace(0, 1, 500)

xmins = np.zeros(len(kappas))
xmaxs = np.zeros(len(kappas))
eps = np.zeros(len(kappas))

for i in range(len(kappas)):
    print(i)

    k = np.load('kappa_'+str(kappas[i])+'.npy')
    sig_data[i,:] = k[120:450]

    
    xmin = np.where(sig_data[i,:] < 0.05)[0]
    xmax = np.where(sig_data[i, :] > 0.95)[0]

    if len(xmin) == 0:
        xmin = np.where(sig_data[i,:] == sig_data[i,:].min())[0][0]
    if len(xmax) == 0:
        xmax = np.where(sig_data[i,:] == sig_data[i,:].max())[0][0]

    else:
        xmin = np.where(sig_data[i,:] < 0.05)[0][-5]
        xmax = np.where(sig_data[i, :] > 0.95)[0][5]

   

    xmins[i] = xmin
    xmaxs[i] = xmax

    eps[i] = x[xmax]-x[xmin]



# fig, ax = plt.subplots()

# for i in range(len(kappas)):

#     if i == 1:
#         ax.plot(sig_data[i,:])
#         ax.axvline(x = xmins[i])
#         ax.axvline(x = xmaxs[i])

# plt.show()

def func(x, a):
    return a * np.sqrt(x)


xdata = np.linspace(0, 1.05*kappas[-1], 100)


popt, pcov = curve_fit(func, kappas, eps)
perr = np.sqrt(np.diag(pcov))

fig = plt.figure()
plt.plot(kappas, eps, 'o', label = r'SpinDec2 Results')
plt.plot(xdata, func(xdata, *popt), 'r-', label = r'$\xi = $'+str(np.round(popt[0], 2))+'$\sqrt{\kappa}$', linestyle = '--')
# plt.fill_between(xdata, func(xdata, *popt) - 2*perr, func(xdata, *popt) + 2*perr, color = 'r', alpha = 0.25)
plt.ylabel(r'Interfacial Width - $\xi$')
plt.xlabel(r'$\kappa$')
# plt.yscale('log')
# plt.xscale('log')
plt.legend(loc  = 'upper left', frameon = False, fontsize = 18)
plt.ylim(0.0, 0.55)
plt.xlim(0, 1.05*kappas[-1])
# plt.fill_between(xdata, func(xdata, *popt) - perr, color = 'r', alpha = 0.25)
plt.show()

print(popt, pcov)



print(perr)
