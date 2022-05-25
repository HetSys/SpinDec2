#!/usr/bin/env python3

import os.path
import numpy as np
from matplotlib import pyplot as plt

def read_timings():
    # Check if timings.txt file exists
    if not os.path.isfile('timings.txt'):
        print('timings.txt file not found')
        exit(1)

    with open('./timings.txt', 'r') as timings:
        times = []
        threads = []

        for line in timings:
            spl = line.split()

            try:
                time = float(spl[0])
            except ValueError:
                continue

            times.append(time)

            if 'omp' == spl[3]:
                pass
                # something to do with associating time as mpi/omp/hybrid
                # Maybe dictionary or another list or something
            elif 'mpi' == spl[3]:
                pass
            elif 'hybrid' in spl[3]:
                pass

            threads.append(spl[-2])


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

def plot_speed_up(problem, speed_up_MPI, speed_up_OMP, speed_up_hybrid, save=false):

    fig1 = plt.figure(figsize=(15,5))

    # MPI
    ax1 = fig1.add_subplot(131)
    ax1.plot(proc_MPI,speed_up_MPI,"-o", label="MPI")
    ax1.set_ylim(1.0,1.4)
    ax1.set_xlabel("no. of processes")

    # OpenMP
    ax2 = fig1.add_subplot(132)
    ax2.plot(proc_OMP,speed_up_OMP,"-o", label="OMP")
    ax2.set_xlabel("no. of processes")

    # Hybrid
    ax3 = fig1.add_subplot(133)
    ax3.plot(MPI_to_OMP,speed_up_hybrid,"-o",label="MPI/OMP")
    ax3.set_ylim(1.0,1.4)
    ax3.set_xlabel("MPI to OMP ratio")

    axes = [ax1, ax2, ax3]

    for ax in axes:
        ax.legend()
        ax.set_ylabel(r"$\psi$")

    problem = problem + r" - Speed up, $\psi$"
    fig1.suptitle(problem)
    
    if save:
        fig1.savefig("speed_up.png", dpi=350, format="png",bbox_inches ="tight")
    
    plt.show()

def plot_prallel_eff(problem, parallel_eff_MPI, parallel_eff_OMP, parallel_eff_hybrid, save=false):

    fig1 = plt.figure(figsize=(15,5))

    # MPI
    ax1 = fig1.add_subplot(131)
    ax1.plot(proc_MPI,parallel_eff_MPI,"-o", c="tab:red", label="MPI")
    ax1.set_ylim(0.0,1.0)
    ax1.set_xlabel("no. of processes")

    # OpenMP
    ax2 = fig1.add_subplot(132)
    ax2.plot(proc_OMP,parallel_eff_OMP,"-o", c="tab:red", label="OMP")
    ax2.set_ylim(0.0,1.0)
    ax2.set_xlabel("no. of processes")

    # Hybrid
    ax3 = fig1.add_subplot(133)
    ax3.plot(proc_hybrid,parallel_eff_hybrid,"-o", c="tab:red", label="MPI/OMP")
    ax3.set_ylim(0.0,1.0)
    ax3.set_xlabel("MPI to OMP ratio")

    axes = [ax1, ax2, ax3]

    for ax in axes:
        ax.legend()
        ax.set_ylabel(r"$\epsilon$")

    problem = problem + r" - Parallel Efficiency, $\epsilon$"
    fig1.suptitle(problem)
    
    if save:
        fig1.savefig("parallel.png", dpi=350, format="png",bbox_inches ="tight")
    
    plt.show()

def plot_karp_flatt(problem, karp_flatt_MPI, karp_flatt_OMP, karp_flatt_hybrid, save=false):

    fig1 = plt.figure(figsize=(15,5))

    # MPI
    ax1 = fig1.add_subplot(131)
    ax1.plot(proc_MPI,karp_flatt_MPI,"-o", c="tab:green", label="MPI")
    ax1.set_ylim(0.0,1.0)
    ax1.set_xlabel("no. of processes")

    # OpenMP
    ax2 = fig1.add_subplot(132)
    ax2.plot(proc_OMP,karp_flatt_OMP,"-o", c="tab:green", label="OMP")
    # ax2.set_ylim(0.0,1.0)
    ax2.set_xlabel("no. of processes")

    # Hybrid
    ax3 = fig1.add_subplot(133)
    ax3.plot(proc_hybrid,karp_flatt_hybrid,"-o", c="tab:green", label="MPI/OMP")
    ax3.set_ylim(0.0,1.0)
    ax3.set_xlabel("MPI to OMP ratio")

    axes = [ax1, ax2, ax3]

    for ax in axes:
        ax.legend()
        ax.set_ylabel("e")

    problem = problem + " - Karp-Flatt Metric, e"
    fig1.suptitle(problem)
    
    if save:
        fig1.savefig("parallel.png", dpi=350, format="png",bbox_inches ="tight")
    
    plt.show()

problem = "Constant Mobilities"

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

# Plot speed up
plot_speed_up(problem,speed_up_MPI, speed_up_OMP, speed_up_hybrid, save=false)

# Plot Parallel Efficiency
plot_prallel_eff(problem, parallel_eff_MPI, parallel_eff_OMP, parallel_eff_hybrid, save=false)

# Plot Karp-Flatt metric
plot_karp_flatt(problem, karp_flatt_MPI, karp_flatt_OMP, karp_flatt_hybrid, save=false)
