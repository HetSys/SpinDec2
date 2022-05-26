#!/usr/bin/env python3

import os.path
import numpy as np
from matplotlib import pyplot as plt

# Set Problem
problem = "Constant Mobilities"

# Timings file to read
input_file = "constant_timings.txt"

def read_timings():
    # Check if timings.txt file exists
    if not os.path.isfile(input_file):
        print(input_file, ' file not found')
        exit(1)

    with open(input_file, 'r') as timings:
        times = []

        for line in timings:
            spl = line.split()

            try:
                time = float(spl[0])
            except ValueError:
                continue

            times.append(time)

    time_hybrid = times[0:4]
    time_mpi = times[4:8]
    serial = times[8]
    time_omp = times[9:]

    return time_hybrid, time_mpi, time_omp, serial


def get_speed_up(serial_time, timings):

    speed_up = np.array([serial_time/timings[i] for i in range(len(timings))])

    return speed_up 

def get_parallel_eff(speed_up,proc):

    parallel_eff = np.array([speed_up[i]/proc[i] for i in range(len(speed_up))])

    return parallel_eff

def get_karp_flatt(speed_up,proc):
    
    e = np.array([((1.0/speed_up[i])-(1.0/proc[i]))/(1.0-(1.0/proc[i])) for i in range(len(proc))])

    return e

def plot_metrics(problem, save=False):

    fig1 = plt.figure(figsize=(20,5))

    # Speed up
    ax1 = fig1.add_subplot(131)
    ax1.plot(proc_MPI,speed_up_MPI,"-o", c="tab:blue", label="MPI")
    ax1.plot(proc_OMP,speed_up_OMP,"-o", c="tab:red", label="OMP")
    ax1.plot(proc_hybrid,speed_up_hybrid,"-o",c="tab:green", label="MPI/OMP")
    ax1.set_ylabel(r"$\psi$")
    ax1.set_ylim(bottom=0.0)

    # Parallel efficiency
    ax2 = fig1.add_subplot(132)
    ax2.plot(proc_MPI,parallel_eff_MPI,"-o", c="tab:blue", label="MPI")
    ax2.plot(proc_OMP,parallel_eff_OMP,"-o", c="tab:red", label="OMP")
    ax2.plot(proc_hybrid,parallel_eff_hybrid,"-o", c="tab:green", label="MPI/OMP")
    ax2.set_ylabel(r"$\epsilon$")
    ax2.set_ylim(0.0,1.0)
   
    # Karp-Flatt metric
    ax3 = fig1.add_subplot(133)
    ax3.plot(proc_MPI,karp_flatt_MPI,"-o", c="tab:blue", label="MPI")
    ax3.plot(proc_OMP,karp_flatt_OMP,"-o", c="tab:red", label="OMP")
    ax3.plot(proc_hybrid,karp_flatt_hybrid,"-o", c="tab:green", label="MPI/OMP")
    ax3.set_ylabel("e")
    ax3.set_ylim(bottom=0.0)

    
    axes = [ax1, ax2, ax3]

    for ax in axes:
        #ax.legend(loc=1, prop={'size': 8}, framealpha=0.2)
        ax.set_xlabel("no. of processes")
        handles, labels = ax.get_legend_handles_labels()
    
    fig1.legend(handles, labels, loc='upper right', prop={"size":12})

    problem = problem + " - Parallel Performance Metrics"
    fig1.suptitle(problem, fontsize=14)
    
    if save:
        fig1.savefig("metrics.png", dpi=350, format="png",bbox_inches ="tight")
    
    plt.show()

MPI_to_OMP = np.array([8, 2, 1, 0.5])
proc_MPI = np.array([4, 16, 25, 36])
proc_OMP = np.array([2, 4, 8, 16, 32, 40])
proc_hybrid = np.array([8, 16, 32, 32])

# Get times
time_hybrid, time_mpi, time_omp, serial = read_timings()

# Rearange times as needed
time_mpi = np.insert(time_mpi, 0, time_mpi[-1])
time_mpi = np.delete(time_mpi, -1)

time_omp = np.insert(time_omp, 0, time_omp[1])
time_omp = np.delete(time_omp, 2)
time_omp = np.insert(time_omp, 1, time_omp[3])
time_omp = np.delete(time_omp, 4)
time_omp = np.insert(time_omp, 2, time_omp[-1])
time_omp = np.delete(time_omp, -1)

time_hybrid = np.insert(time_hybrid,4,time_hybrid[0])
time_hybrid = np.delete(time_hybrid, 0)

# Get speed up
speed_up_MPI = get_speed_up(serial,time_mpi)
speed_up_OMP = get_speed_up(serial,time_omp)
speed_up_hybrid = get_speed_up(serial,time_hybrid)

# Get parallel efficiency
parallel_eff_MPI = get_parallel_eff(speed_up_MPI, proc_MPI)
parallel_eff_OMP = get_parallel_eff(speed_up_OMP, proc_OMP)
parallel_eff_hybrid = get_parallel_eff(speed_up_hybrid, proc_hybrid)

# Get karp-Flatt Metric
karp_flatt_MPI = get_karp_flatt(speed_up_MPI, proc_MPI)
karp_flatt_OMP = get_karp_flatt(speed_up_OMP, proc_OMP)
karp_flatt_hybrid = get_karp_flatt(speed_up_hybrid, proc_hybrid)

# Plot metrics
plot_metrics(problem, save=True)
