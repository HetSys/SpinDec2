import netCDF4 as NC
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from scipy.ndimage import uniform_filter1d
import scipy.optimize as cf
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse


def f(x,b,a):
    return b*x+a


def nbs(p,c):
    ns = []
    s = c.shape
    for i in [-1,0,1]:
        for j in [-1,0,1]:
            if(i !=0 or j !=0):
                px = p[0]+i
                py = p[1]+j
                if(px > s[0]-1):
                    px -= s[0]
                if(px < 0):
                    px+=s[0]
                if(py > s[0]-1):
                    py -= s[0]
                if(py < 0):
                    py+=s[0]
                if c[px,py] < 1e-2:
                    ns.append([px,py])
    return ns

def get_area(sp,c,counted):
    if c[sp[0],sp[1]] >1e-2:
        return 1
    to_check = [sp]
    checked = []
    counted.append(sp)
    area = 1
    while to_check:
        checked.append(to_check[0])
        for i in nbs(to_check[0],c):
            if i not in checked and i not in counted:
                to_check.append(i)
            if i not in counted:
                area += 1
                #c[i[0],i[1]]+=1
                counted.append(i)
        del to_check[0]
    return area




def get_domain_szie(c):
    lens = 0
    tot_lens = 0
    no_lens = 0
    for j in c:
        for i in j:
            if i>0.7:
                lens = lens +1
            else:
                no_lens = no_lens +1
                tot_lens = tot_lens + lens
                lens = 0
        for i in j:
            if i<0.3:
                lens = lens +1
            else:
                no_lens = no_lens +1
                tot_lens = tot_lens + lens
                lens = 0
    return tot_lens/no_lens

data = NC.Dataset("CH_output.nc", "r", format="NETCDF4")

c = np.array([data.variables['c'][:]][0])


#95,33

sizes = []
for i in c:
    sizes.append(get_domain_szie(i))

x_dat = np.linspace(0,0.1,1000)


mens = []
for t in range(300):
    areas = []
    diffc = np.power(np.diff(c[t],axis=0),2)[:,:-1]+np.power(np.diff(c[t],axis=1)[:-1,:],2)
    #plt.imshow(diffc)
    #plt.show()
    counted = []
    #diffc[33,95]+=1
    #for i in range(diffc.shape[0]):
        #for j in range(diffc.shape[0]):
    area = get_area([26,44],diffc,counted)
            #if area != 1:
    areas.append(area)
    print(np.mean(areas),t)
    mens.append(np.mean(areas))
s = c.shape
mens = np.array(mens)
cofs , _ = cf.curve_fit(f,np.log(x_dat[50:300]),np.log(mens[50:300]))
print(cofs)
plt.plot(np.log(x_dat[:300]),f(np.log(x_dat[:300]),*cofs),label = "fit")
plt.plot(np.log(x_dat[:300]),np.log(mens),label="Size of domain containing point p")
plt.legend()
plt.xlabel("log(t)")
plt.ylabel("log(domain size)")
plt.show()
plt.plot(mens/(s[0]*s[1]))
plt.xlabel("t(s)")
plt.ylabel("domain size")
plt.show()
