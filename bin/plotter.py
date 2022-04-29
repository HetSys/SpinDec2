#!/usr/bin/env python3

import netCDF4 as NC
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from scipy.ndimage import uniform_filter1d
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm


# Pre-set plotting style 

plt.rcParams["figure.figsize"] = [8, 8]
plt.rcParams["figure.autolayout"] = True
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['legend.fontsize'] = 14


# simple moving average using scipy uniform_filter1d pacakge
# smooths out free energy result
def moving_av(Q, span, period = 10):
    return uniform_filter1d(Q, size = span // period)

#calcualtes the value of f(c = c_av) (average bulk free energy)
def bulk_energy_av(coeffs, c_av):
    bulk = 0
    for i in range(len(coeffs)):
        bulk += coeffs[i] * c_av**i
    return bulk

#Visualisation Code Begins
class Vis_CH:

    """Class implementing Visualisation of Cahn-Hilliard"""

    #initilisation of the class
    def __init__(self, file, av_period):

        self.file = file
        self.av_period = av_period

    #read in data from output file (netcdf)
    def read_netcdf(self, verbose = True):

        data = NC.Dataset(self.file, "r", format="NETCDF4")
    
        #concentration order parameter c(t, y, x)
        self.c = [data.variables['c'][:]][0]

        #### Issue/needs checking - does c need to be transposed?
        
        #free energy over time
        self.F_tot = [data.variables['F_tot'][:]][0]

        #user inputted coefficients for f(c)
        self.coeffs = [data.variables['coeffs'][:]][0]

        
        self.dt = data.dt         #time-step
        self.N_t = int(data.N_t)  #time iterations
        self.N_x = int(data.N_x)  #x dicretizations
        self.N_y = int(data.N_y)  #y discretizations
        self.M_A = data.M_A       #atomic mobility (species A)
        self.M_B = data.M_B       #atomic mobility (species B)
        self.c_0 = data.c_0       #user set initail grid average c
        self.kappa = data.kappa   #gradient term coefficient

        #time array
        self.time = np.arange(self.N_t) * self.dt

        #moving average of F
        if (self.av_period > self.N_t):
            print("WARNING:")
            print('set average period is greater than total time period (N_t = ' + str(self.N_t) + ' setting to default value of 10')
            self.av_period = 10

        self.F_av = moving_av(self.F_tot, self.N_t, self.av_period)

        #frame speed of animation depending on number of files (i.e more files, want to be quicker)

        if (self.N_t >= 1000):
            self.frame_step = 1
        else:
            self.frame_step = 5

        #print metadata to terminal
        if (verbose):

            print('\u0394t =', self.dt, 's')
            print('N_x, N_y, N_t =', self.N_x, self.N_y, self.N_t )
            print('\u039C_A =', self.M_A)
            print('\u039C_B =', self.M_B)
            print('\u03BA =', self.kappa)
            print('c_0', self.c_0)
            string = 'f(c) = ' + str(np.round(self.coeffs[0], 2)) + ' + '
            for i in range(len(self.coeffs)):
                if i == 0:
                    continue
                if i == len(self.coeffs)-1:
                    string += str(np.round(self.coeffs[i],2)) + 'c^'+str(i)
                    continue
                if (self.coeffs[i+1] < 0):
                    string += str(np.round(self.coeffs[i],2)) + 'c^'+str(i) + ' '
                else:
                    string += str(np.round(self.coeffs[i],2)) + 'c^'+str(i) + ' + '

            print(string)

    #function that plots an animation of both c and F (and its moving average) on the same figure
    def dual_animation(self):
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15,15))
        cmap = cm.jet
        im = ax1.imshow(self.c[0], interpolation='none', aspect='auto', vmin=0, vmax=1, cmap = cmap)
        fig.suptitle(r"Time =  " + str(self.time[0]))
        div = make_axes_locatable(ax1)
        cax = div.append_axes('right', '5%', '5%')
        fig.colorbar(im, cax = cax)

        def animate_grid(i):

            im.set_array(self.c[i])
            fig.suptitle(r"Time =  " + str(self.time[i]))
            ax1.set_xlabel(r'x')
            ax1.set_ylabel(r"$y$")
            ax1.set_title(r'$c(x,y)$')
            div = make_axes_locatable(ax1)
            cax = div.append_axes('right', '5%', '5%')
            fig.colorbar(im, cax = cax)

        

        self.x = [self.time[0]]
        self.y = [self.F_tot[0]]
        self.y_av = [self.F_av[0]]

        ax2.plot(self.x,self.y, scaley=True, scalex=True, color="blue", label = r'CH Data')
        ax2.plot(self.x,self.y_av, scaley=True, scalex=True, color="red", linestyle = '--', label = r'Moving Average ($t_{period}$ = ' + str(self.N_t // self.av_period) + r' s)')
        ax2.set_xlabel(r'Time (s)')
        ax2.set_ylabel(r"Total Free Energy")
        ax2.set_title(r'$F(t)$')
        ax2.legend(loc = 'upper left')

        def animate_line(j):
            j += 1
            self.x.append(self.time[j])
            self.y.append(self.F_tot[j])
            self.y_av.append(self.F_av[j])

            ax2.plot(self.x,self.y, scaley=True, scalex=True, color="blue")
            ax2.plot(self.x,self.y_av, scaley=True, scalex=True, color="red", linestyle = '--')

        anim_grid = animation.FuncAnimation(fig = fig, func = animate_grid, frames=len(self.time), repeat=False, interval=self.frame_step)
        anim_line = animation.FuncAnimation(fig = fig, func = animate_line, frames = len(self.time) - 1, repeat = False, interval=self.frame_step)
        plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.5, 
                        hspace=0.4)
        plt.show()
 
    #function that plots an animation of c
    def grid_animation(self):
        
        fig, ax = plt.subplots()
        cmap = cm.jet
        im = ax.imshow(self.c[0], interpolation='none', aspect='auto', vmin=0, vmax=1, cmap = cmap)
        fig.suptitle(r"Time =  " + str(self.time[0]))
        div = make_axes_locatable(ax)
        cax = div.append_axes('right', '5%', '5%')
        fig.colorbar(im, cax = cax)

        def animate_grid(i):

            im.set_array(self.c[i])
            fig.suptitle(r"Time =  " + str(self.time[i]))
            ax.set_xlabel(r'x')
            ax.set_ylabel(r"$y$")
            ax.set_title(r'$c(x,y)$')
            div = make_axes_locatable(ax)
            cax = div.append_axes('right', '5%', '5%')
            fig.colorbar(im, cax = cax)

        anim_grid = animation.FuncAnimation(fig = fig, func = animate_grid, frames=len(self.time), repeat=False, interval=self.frame_step)
        plt.show()

    #function that plots an animation of F (and its moving average)
    def F_animation(self):
        fig, ax = plt.subplots()

        self.x = [self.time[0]]
        self.y = [self.F_tot[0]]
        self.y_av = [self.F_av[0]]

        ax.plot(self.x,self.y, scaley=True, scalex=True, color="blue", label = r'CH Data')
        ax.plot(self.x,self.y_av, scaley=True, scalex=True, color="red", linestyle = '--', label = r'Moving Average ($t_{period}$ = ' + str(self.N_t // self.av_period) + r' s)')
        ax.set_xlabel(r'Time (s)')
        ax.set_ylabel(r"Total Free Energy")
        ax.set_title(r'$F(t)$')
        ax.legend(loc = 'upper left')

        def animate_line(j):
            j += 1
            self.x.append(self.time[j])
            self.y.append(self.F_tot[j])
            self.y_av.append(self.F_av[j])

            ax.plot(self.x,self.y, scaley=True, scalex=True, color="blue")
            ax.plot(self.x,self.y_av, scaley=True, scalex=True, color="red", linestyle = '--')

        anim_line = animation.FuncAnimation(fig = fig, func = animate_line, frames = len(self.time) - 1, repeat = False, interval=self.frame_step)
        plt.show()

    #Function that plots F (and its moving average)
    def F_plot(self):
        fig, ax = plt.subplots()
        ax.plot(self.time,self.F_tot, scaley=True, scalex=True, color="blue", label = r'CH Data')
        ax.plot(self.time,self.F_av, scaley=True, scalex=True, color="red", linestyle = '--', label = r'Moving Average ($t_{period}$ = ' + str(self.N_t // self.av_period) + r' s)')
        ax.set_xlabel(r'Time (s)')
        ax.set_ylabel(r"Total Free Energy - $F(t)$")
        ax.legend(loc = 'upper left')
        plt.show()

    #function thats plots the initial bulk free energy and the trajectory of the 
    #averge bulk free energy over time in c c-f(c) space
    def bulk_energy_traj(self):

            self.c_av = np.array([])
            self.f_b_av = np.array([])
            self.c_space = np.linspace(-0.5, 1.5, 100)
            self.f_b = np.zeros(len(self.c_space))

            for i in range(len(self.coeffs)):
                self.f_b += self.coeffs[i] * self.c_space**i

            for i in range(self.N_t):
                av = np.average(self.c[i])
                self.c_av = np.append(self.c_av, av)
                self.f_b_av = np.append(self.f_b_av, bulk_energy_av(self.coeffs, self.c_av))

            self.traj_c = [self.c_av[0]]
            self.traj_f_b = [self.f_b_av[0]]

            fig, ax = plt.subplots()

            ax.plot(self.c_space, self.f_b, label = r'$f(c)$', color = 'blue', linestyle = '--')
            fig.suptitle(r"Time =  " + str(self.time[0]))
            ax.scatter(self.traj_c, self.traj_f_b, label = r'Trajectory', color = 'red')
            ax.set_xlabel(r'c')
            ax.set_ylabel(r'Bulk Free Energy')
        
            def animate_traj(i):
                i += 1
                self.traj_c.append(self.c_av[i])
                self.traj_f_b.append(self.f_b_av[i])
                ax.plot(self.c_space, self.f_b, color = 'blue', linestyle = '--')
                fig.suptitle(r"Time =  " + str(self.time[i]))
                ax.scatter(self.traj_c, self.traj_f_b, color = 'red')

            anim_traj = animation.FuncAnimation(fig = fig, func = animate_traj, frames = len(self.time) - 1, repeat = False, interval=self.frame_step)
        
            plt.show()
    

#run visualisation

#perhaps need to change to get the file from whatever dir we choose
res_CH = Vis_CH(file = 'CH_output.nc', av_period=10)

res_CH.read_netcdf(verbose=True)

res_CH.dual_animation()

res_CH.grid_animation()

res_CH.F_animation()

res_CH.F_plot()

res_CH.bulk_energy_traj()
