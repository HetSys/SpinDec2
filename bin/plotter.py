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
    def read_and_setup(self, verbose = True):

        data = NC.Dataset(self.file, "r", format="NETCDF4")

        #concentration order parameter c(t, y, x)
        self.c = np.array([data.variables['c'][:]][0])
        # self.c = np.swapaxes(self.c, 1, 2)

        #### Issue/needs checking - does c need to be transposed?

        #free energy over time
        self.F_tot = np.array([data.variables['F_tot'][:]][0])

        #user inputted coefficients for f(c)
        self.coeffs = np.array([data.variables['coeffs'][:]][0])


        self.dt = data.dt         #time-step
        self.N_t = int(data.Nt)  #time iterations
        self.N_x = int(data.Nx)  #x dicretizations
        self.N_y = int(data.Ny)  #y discretizations
        self.M_A = data.MA       #atomic mobility (species A)
        self.M_B = data.MB       #atomic mobility (species B)
        self.c_0 = data.c0       #user set initail grid average c
        self.kappa = data.kappa   #gradient term coefficient

        #time array
        self.time = np.arange(self.N_t) * self.dt

         #moving average of F
        if (self.av_period > self.N_t):
            print("WARNING:")
            print('set average period is greater than total time period (N_t = ' + str(self.N_t) + ' setting to default value')

            self.av_period = 10

        self.F_av = moving_av(self.F_tot, self.N_t, self.av_period)

        #print metadata to terminal
        if (verbose):

            print('\u0394t =', self.dt, 's')
            print('N_x, N_y, N_t =', self.N_x, self.N_y, self.N_t )
            print('\u039C_A =', self.M_A)
            print('\u039C_B =', self.M_B)
            print('\u03BA =', self.kappa)
            print('c_0 =', self.c_0)
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

        

        self.nframes = 300

        if self.N_t < self.nframes:
            self.modulo = 1
        else:
            self.modulo = int(self.N_t/self.nframes)

        c_plot = []
        self.F_anim = np.array([])
        self.F_anim_av = np.array([])
        self.time_anim = np.array([])

        for i in range(len(self.c)):
            if (i % self.modulo == 0):
                c_plot.append(self.c[i])
                self.F_anim = np.append(self.F_anim, self.F_tot[i])
                self.F_anim_av = np.append(self.F_anim_av, self.F_av[i])
                self.time_anim = np.append(self.time_anim, self.time[i])

        self.c_plot = np.array(c_plot)

        
        

    
    #function that plots an animation of both c and F (and its moving average) on the same figure
    def dual_animation(self):

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15,15))

        # cmap = cm.jet
        pos1 = ax1.imshow(np.zeros(shape=(self.N_y, self.N_x)), vmin=0, vmax=1, origin='lower')
        div = make_axes_locatable(ax1)
        cax = div.append_axes('right', '5%', '5%')
        fig.colorbar(pos1, cax = cax)

        ax1.set_xlabel(r'x')
        ax1.set_ylabel(r"$y$")
        ax1.set_title(r'$c(x,y)$')

        line, = ax2.plot([], [], '-', label = r'Sim Data', color = 'black')
        line_av, = ax2.plot([], [], '--', label = r'Moving Average ($t_{period}$ = ' + str(self.N_t // self.av_period) + r' s)', color = 'red')
        ax2.set_xlabel(r'Time (ns)')
        ax2.set_ylabel(r'Total Free Energy - $F$')
        ax2.set_xlim(self.time_anim[0], self.time_anim[-1])
        ax2.set_ylim(0.98*self.F_anim[0], 1.02*self.F_anim[-1])
        ax2.legend()

        def init():
            pos1.set_array(self.c_plot[0])
            line.set_data(self.time_anim[0], self.F_anim[0])
            line_av.set_data(self.time_anim[0], self.F_anim_av[0])
            return [pos1, line, line_av,]

        def animate(i):
            pos1.set_array(self.c_plot[i])
            line.set_data(self.time_anim[:i], self.F_anim[:i])
            line_av.set_data(self.time_anim[:i], self.F_anim_av[:i])
            return [pos1, line, line_av,]
    
        anim = animation.FuncAnimation(fig = fig, func = animate, init_func = init, frames = len(self.c_plot), repeat = False, interval=1,blit=True)

        plt.show()

       


    #function that plots an animation of both c and F (and its moving average) on the same figure
    def F_animation(self):

        fig, ax = plt.subplots()
        
        line, = ax.plot([], [], '-', label = r'Sim Data', color = 'black')
        line_av, = ax.plot([], [], '--', label = r'Moving Average ($t_{period}$ = ' + str(self.N_t // self.av_period) + r' s)', color = 'red')
        ax.set_xlabel(r'Time (ns)')
        ax.set_ylabel(r'Total Free Energy - $F$')
        ax.set_xlim(self.time_anim[0], self.time_anim[-1])
        ax.set_ylim(0.98*self.F_anim[0], 1.02*self.F_anim[-1])
        ax.legend()

        def init():
            line.set_data(self.time_anim[0], self.F_anim[0])
            line_av.set_data(self.time_anim[0], self.F_anim_av[0])
            return [line,line_av,]

        def animate(i):
            line.set_data(self.time_anim[:i], self.F_anim[:i])
            line_av.set_data(self.time_anim[:i], self.F_anim_av[:i])
            return [line,line_av,]

        anim_traj = animation.FuncAnimation(fig = fig, func = animate, init_func = init, frames = len(self.c_plot), repeat = False, interval=1,blit=True)
               

        plt.show()
    
    #function that plots an animation of c
    def grid_animation(self):

        fig, ax = plt.subplots()
        # cmap = cm.jet
           # cmap = cm.jet
        pos1 = ax.imshow(np.zeros(shape=(self.N_y, self.N_x)), vmin=0, vmax=1, origin='lower')
        div = make_axes_locatable(ax)
        cax = div.append_axes('right', '5%', '5%')
        fig.colorbar(pos1, cax = cax)

        ax.set_xlabel(r'x')
        ax.set_ylabel(r"$y$")
        ax.set_title(r'$c(x,y)$')

        def init():
            pos1.set_array(self.c_plot[0])
            return [pos1]

        def animate(i):
            pos1.set_array(self.c_plot[i])
            return [pos1]
    
        anim = animation.FuncAnimation(fig = fig, func = animate, init_func = init, frames = len(self.c_plot), repeat = False, interval=1,blit=True)

        plt.show()


    #function that plots an animation of c
    def grid_snapshot(self, snapshot = -1):

        fig, ax = plt.subplots()
        # cmap = cm.jet
        pos1 = ax.imshow(self.c_plot[snapshot][:,:], vmin=0, vmax=1, origin = 'lower')
        div = make_axes_locatable(ax)
        cax = div.append_axes('right', '5%', '5%')
        fig.colorbar(pos1, cax = cax)

        ax.set_xlabel(r'x')
        ax.set_ylabel(r"$y$")
        ax.set_title(r'$c(x,y)$')

        plt.show()


    #Function that plots F (and its moving average)
    def F_plot(self):
        fig, ax = plt.subplots()
        ax.plot(self.time,self.F_tot, scaley=True, scalex=True, color="blue", label = r'CH Data')
        ax.plot(self.time,self.F_av, scaley=True, scalex=True, color="red", label = r'Moving Average ($t_{period}$ = ' + str(self.N_t // self.av_period) + r' s)')
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

            for i in range(len(self.c_plot)):
                av = np.average(self.c_plot[i])
                self.c_av = np.append(self.c_av, av)
                self.f_b_av = np.append(self.f_b_av, bulk_energy_av(self.coeffs, self.c_av))

            fig, ax = plt.subplots()

            ax.plot(self.c_space, self.f_b, label = r'$f(c)$', color = 'blue', linestyle = '--')
            traj, = ax.plot([], [], 'o', label = r'Trajectory', color = 'red')
            ax.set_xlabel(r'c')
            ax.set_ylabel(r'Bulk Free Energy')
            ax.set_xlim(0.49, 0.51)

            def init():
                traj.set_data(self.c_av[0], self.f_b_av[0])
                return traj,

            def animate(i):
                traj.set_data(self.c_av[:i], self.f_b_av[:i])
                return traj,
               
                

            anim_traj = animation.FuncAnimation(fig = fig, func = animate, init_func = init, frames = len(self.c_plot), repeat = False, interval=1,blit=True)

            plt.show()


#run visualisation

#perhaps need to change to get the file from whatever dir we choose
res_CH = Vis_CH(file = 'CH_output.nc', av_period=10)

res_CH.read_and_setup(verbose=True)

# res_CH.F_animation()

res_CH.grid_animation()

res_CH.grid_snapshot(snapshot=-1)

# res_CH.F_plot()

# res_CH.bulk_energy_traj()

# res_CH.dual_animation()