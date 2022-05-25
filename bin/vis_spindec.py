
  
#!/usr/bin/env python3

import netCDF4 as NC
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from scipy.ndimage import uniform_filter1d
from mpl_toolkits.axes_grid1 import make_axes_locatable
import argparse


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

#calcualtes the value of f(c)
def bulk_energy(c, coeffs):
    bulk = 0
    for i in range(len(coeffs)):
        bulk += coeffs[i] * c**i
    return bulk

#calcualtes the gradient of f(c)
def mu(c, coeffs):
    bulk = 0
    for i in range(len(coeffs)):
        bulk += i*coeffs[i] * c**(i-1)
    return bulk

#calcualtes the curvature of f(c)
def dmu_dc(c, coeffs):
    bulk = 0
    for i in range(len(coeffs)):
        bulk += i*(i-1)*coeffs[i] * c**(i-2)
    return bulk


# estimates the positions of the extrema (min and maxima) of f(c)
def bulk_energy_extrema(coeffs):

    c = np.linspace(-2,2,1000000)

    min_pos = np.array([])
    max_pos = np.array([])

    df_dc = np.zeros(len(c))
        
    for i in range(len(c)):

        df_dc = mu(c[i], coeffs)
        d2f_dc2 = dmu_dc(c[i], coeffs)

        if (np.abs(df_dc) < 1e-3 and d2f_dc2 >0):

            min_pos = np.append(min_pos, np.round(c[i],2))

        if (np.abs(df_dc) < 1e-3 and d2f_dc2 <0):

            max_pos = np.append(max_pos, np.round(c[i],2))

    return min_pos, max_pos

# Argument parser for saving and viewing the  plot
parser = argparse.ArgumentParser(description='Visulise spinodal decomposistion in SpinDec2 (DEFAULT: Saves plots)')

# controls snapshot that you can view in grid_snapshot plotter function
parser.add_argument('-sn', '--snap', type=int, default=-1, metavar = '',help='Select Snapshot Of Grid To Plot')

# controls whether user wants to save plots only
parser.add_argument('-s', '--save', action='store_true', help='Save plots only')
# controls whether user wants to view plots only
parser.add_argument('-d', '--display', action='store_true', help='View plots only')


args = parser.parse_args()


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
            print('###### SpinDec2 Meta Data ######')

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

        
        # Number of frames for the aniamtion
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

        # For the bounds of the plots
        self.minima, self.maxima  = bulk_energy_extrema(self.coeffs)


        

    
    #function that plots an animation of both c and F (and its moving average) on the same figure
    def dual_animation(self, save_plot = True):

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15,15))

        vmin = np.min(self.minima)
        vmax = np.max(self.minima)
        pos1 = ax1.imshow(np.zeros(shape=(self.N_y, self.N_x)), vmin=vmin, vmax=vmax, origin='lower')
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
            if (save_plot):
                fig.suptitle('Time = '+str(self.time_anim[0]), fontsize = 20)
            return [pos1, line, line_av,]

        def animate(i):
            pos1.set_array(self.c_plot[i])
            line.set_data(self.time_anim[:i], self.F_anim[:i])
            line_av.set_data(self.time_anim[:i], self.F_anim_av[:i])
            if (save_plot):
                fig.suptitle('Time = '+str(self.time_anim[i]), fontsize = 20)
            return [pos1, line, line_av,]
    
        anim = animation.FuncAnimation(fig = fig, func = animate, init_func = init, frames = len(self.c_plot), repeat = False, interval=1,blit=True)

        if (save_plot):

            FFwriter = animation.FFMpegWriter(fps = 50)
            anim.save('Dual_Animation.mp4', writer = FFwriter)

        else:
            plt.show()

       


    #function that plots an animation of both c and F (and its moving average) on the same figure
    def F_animation(self, save_plot = True):

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

        if (save_plot):
        
            FFwriter = animation.FFMpegWriter(fps = 50)
            anim_traj.save('Free_Energy_Animation.mp4', writer = FFwriter)
               
        else:
            plt.show()
    
    #function that plots an animation of c
    def grid_animation(self, save_plot = True):

        fig, ax = plt.subplots()
        vmin = np.min(self.minima)
        vmax = np.max(self.minima)
        pos1 = ax.imshow(np.zeros(shape=(self.N_y, self.N_x)), vmin=vmin, vmax=vmax, origin='lower')
        div = make_axes_locatable(ax)
        cax = div.append_axes('right', '5%', '5%')
        fig.colorbar(pos1, cax = cax)

        ax.set_xlabel(r'x')
        ax.set_ylabel(r"$y$")
        ax.set_title(r'$c(x,y)$')

        def init():
            pos1.set_array(self.c_plot[0])
            if (save_plot):
                fig.suptitle('Time = '+str(self.time_anim[0]), fontsize = 20)
            return [pos1]

        def animate(i):
            pos1.set_array(self.c_plot[i])
            if (save_plot):
                fig.suptitle('Time = '+str(self.time_anim[i]), fontsize = 20)
            return [pos1]
    
        anim = animation.FuncAnimation(fig = fig, func = animate, init_func = init, frames = len(self.c_plot), repeat = False, interval=1,blit=True)

        if (save_plot):

            FFwriter = animation.FFMpegWriter(fps = 50)
            anim.save('Grid_Animation.mp4', writer = FFwriter)

        else:
            plt.show()


    #function that plots an animation of c
    def grid_snapshot(self, snapshot = -1, save_plot = True):

        fig, ax = plt.subplots()
        vmin = np.min(self.minima)
        vmax = np.max(self.minima)
        pos1 = ax.imshow(self.c_plot[snapshot][:,:], vmin=vmin, vmax=vmax, origin = 'lower')
        div = make_axes_locatable(ax)
        cax = div.append_axes('right', '5%', '5%')
        fig.colorbar(pos1, cax = cax)
        fig.suptitle('Time = '+str(self.time[snapshot]))

        ax.set_xlabel(r'x')
        ax.set_ylabel(r"$y$")
        ax.set_title(r'$c(x,y)$')

        if (save_plot):
            plt.savefig('Grid_Snaphot.png', bbox_inches='tight')
        else:
            plt.show()


    #Function that plots F (and its moving average)
    def F_plot(self, save_plot = True):
        fig, ax = plt.subplots()
        ax.plot(self.time,self.F_tot, scaley=True, scalex=True, color="blue", label = r'CH Data')
        ax.plot(self.time,self.F_av, scaley=True, scalex=True, color="red", label = r'Moving Average ($t_{period}$ = ' + str(self.N_t // self.av_period) + r' s)')
        ax.set_xlabel(r'Time (s)')
        ax.set_ylabel(r"Total Free Energy - $F(t)$")
        ax.legend(loc = 'upper left')

        if (save_plot):
            plt.savefig('Free_Energy.png', bbox_inches='tight')
        else:
            plt.show()

   
    #function thats plots the initial bulk free energy and the trajectory of the
    #averge bulk free energy over time in c c-f(c) space
    def bulk_energy_traj(self, save_plot = True):
        
        if (self.N_x == self.N_y):
            if (self.N_x == 64 or self.N_x == 128 or self.N_x == 256):

            
                self.c_av = np.zeros((len(self.c_plot), 16))
                self.f_b_av = np.zeros((len(self.c_plot), 16))

                self.c_space = np.linspace(-2, 2, 100)
                self.f_b = np.zeros(len(self.c_space))

                for i in range(len(self.coeffs)):
                    self.f_b += self.coeffs[i] * self.c_space**i


                for i in range(len(self.c_plot)):
                    c_split_x = np.array_split(self.c_plot,4, axis=0)
                    c_grid = []
                    for grid in c_split_x :
                        c_split_y = np.array_split(grid,4,axis=1)
                        c_grid += c_split_y
                
                    for j in range(16):
                        av = np.average(c_grid[j])
                        self.c_av[i,j] = av
                        self.f_b_av[i,j] = bulk_energy(av, self.coeffs)

                
                fig, axs = plt.subplots(4, 4, figsize=(20,20))

                y_maxima = np.array([]) 
                y_minima = np.array([])

                for i in range(len(self.minima)):
                    y_minima = np.append(y_minima, bulk_energy(self.minima[i], self.coeffs))

                for i in range(len(self.maxima)):
                    y_maxima = np.append(y_maxima, bulk_energy(self.maxima[i], self.coeffs))

          
                count = 0
                for i in range(4):
                    for j in range(4):
                        axs[i,j].plot(self.c_space, self.f_b, label = r'$f(c)$', color = 'blue', linestyle = '--')
                        axs[i,j].plot(self.c_av[:,count], self.f_b_av[:,count], 'o', label = 'Traj', color = 'r')

                        count += 1

                        axs[i,j].set_xlabel(r'c')
                        axs[i,j].set_ylabel(r'Bulk Free Energy')
                        axs[i,j].legend()
                        if self.minima.min() == 0 or self.minima.max() == 0:
                            axs[i,j].set_xlim(-0.5, 1.5)
                        else:
                            axs[i,j].set_xlim(0.85*self.minima.min(), 1.15*self.minima.max())

                        axs[i,j].set_ylim(0.95*y_minima.min(), 1.20*y_maxima.max())

                
                if (save_plot):
                    plt.savefig('Bulk_Energy_Traj.png', bbox_inches='tight')
                else:
                    plt.show()

            
            else:
                print('Function Terminated, Error = Size: Plotting Routine (bulk_energy_traj) can only be used for sqaure grids of size (64,64) ; (128, 128) or (256, 256)')
                return None

        else:
            print('Function Terminated, Error = Shape: Plotting Routine (bulk_energy_traj) can only be used for sqaure grids of size (64,64) ; (128, 128) or (256, 256)')
            return None
                

#run visualisation

#perhaps need to change to get the file from whatever dir we choose
res_CH = Vis_CH(file = 'CH_output.nc', av_period=10)

res_CH.read_and_setup(verbose=True)

# arg parser for saving or viewing
if args.save:
    save = True
if args.display:
    save = False
else:
    save = True

res_CH.F_animation(save_plot=save)
if save:
    print('Saved Free Energy Animation')

res_CH.grid_animation(save_plot=save)
if save:
    print('Saved Grid Animation')

res_CH.grid_snapshot(snapshot=args.snap, save_plot=save)
if save:
    print('Saved Grid Snapshot N_t =' + str(args.snap))

res_CH.F_plot(save_plot=save)
if save:
    print('Saved Free Energy plot')

res_CH.bulk_energy_traj(save_plot=save)
if save:
    print('Saved Bulk Energy Trajectory Plot')

res_CH.dual_animation(save_plot=save)
if save:
    print('Saved Dual Animation')