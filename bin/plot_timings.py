#!/usr/bin/env python3

import os.path


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


def plot_timings():
    pass


if __name__ == '__main__':
    read_timings()
    plot_timings()
