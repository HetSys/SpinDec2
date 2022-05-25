#!/bin/bash
# Automate testing in parallel for SpinDec
# Usage: ./auto_profile <dir_name> <problem_type>

trap "exit" INT

prepare_profile () {
    ### Create dirs to run spindec ###
    cd "$1" || exit 1

    threads=("$@")
    threads=("${threads[@]/$1}")

    for i in "${threads[@]:1}"; do
        mkdir "$i"
        cd "$i" || exit 1
        cp ../../input.txt ./
        cd ../
    done

    cd ../
}

# Check if user supplied directory to perform computations in
if [[ -z $1 || -z $2 ]]; then
    echo "Required argument(s) not specified"
    exit 1
fi

# Threads to run MPI and OpenMP with
omp_threads=(1 2 4 8 16 32 40)
mpi_threads=(4 16 25 36)
hybrid_mpi_4_threads=(2 4 8)
hybrid_mpi_16_threads=(2)

# Load required modules
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 netCDF-Fortran/4.4.4 FFTW/3.3.8

mkdir "$1" && cd "$1" || exit 1

# Create input file to use
echo "Concentration_max = 0.9
Concentration_min = 0.1
Domain_x_size = 600
Domain_y_size = 600
Mobility_A = 4
Mobility_B = 4
free_energy_gradient_parameter = 0.0001
Bulk_free_energy = 1
Checkpointing_interval = 5000000
Max_time = 1e-5
time_step = 1e-8
dF_tolerance = 2.0
Random_seed = 12345356
Use_input = 0
Exitation_A = 0.1
Exitation_B = 0.2
Temperature_min = 800
Temperature_max = 1000
Problem = "$2"
Stabilization_Term = -1" > input.txt

# Create test dirs
mkdir omp mpi hybrid_mpi_4 hybrid_mpi_16

# Create dirs for each profile case
prepare_profile omp "${omp_threads[@]}"
prepare_profile mpi "${mpi_threads[@]}"
prepare_profile hybrid_mpi_4 "${hybrid_mpi_4_threads[@]}"
prepare_profile hybrid_mpi_16 "${hybrid_mpi_16_threads[@]}"

rm input.txt

for i in *; do (
    cd "$i" || exit 1
    for j in *; do (
        cd "$j" || exit 1

        # Run the computations and print timings to timings.txt
        if [[ "$i" == *"omp"* ]]; then
            export OMP_NUM_THREADS="$j" && mpirun -np 1 spindec | tee prof_out.txt
        elif [[ "$i" == "mpi" ]]; then
            export OMP_NUM_THREADS=1 && mpirun -np "$j" spindec | tee prof_out.txt
        elif [[ "$i" == *"mpi_4"* ]]; then
            export OMP_NUM_THREADS="$j" && mpirun -np 4 spindec | tee prof_out.txt
        elif [[ "$i" == *"mpi_16"* ]]; then
            export OMP_NUM_THREADS="$j" && mpirun -np 16 spindec | tee prof_out.txt
        else
            echo "An error occurred"
            exit 1
        fi

        # Parse timings for each run
        # Timings are written to ./$1/timings.txt
        # If anyone sees this, I'm pretty proud of this command
        awk -v type="$i" -v num="$j" -F ' ' 'END{ print $3, "sec on " type " with " num " threads" }' \
            prof_out.txt >> ../../timings.txt

    ) done
) done
