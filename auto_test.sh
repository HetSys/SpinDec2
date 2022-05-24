#!/bin/bash

omp_threads=(2 4 8 16 32 40)
mpi_threads=(4 16 25 36)
hybrid_mpi_4_threads=(2 4 8)
hybrid_mpi_16_threads=(2)

if [[ -z $1 ]]; then
    echo "Directory argument not specified"
    exit 1
fi

mkdir $1 && cd $1
echo 'Concentration_max = 0.9
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
Problem = "Constant"
Stabilization_Term = -1' > input.txt

mkdir omp mpi hybrid_mpi_4 hybrid_mpi_16

prepare_profile () {
    cd "$1" || exit 1

    for i in "$2"; do
        mkdir "$i"
        ( cd "$i" && cp ../../input.txt ./)
    done
}

prepare_profile omp "${omp_threads[@]}"
prepare_profile mpi "${mpi_threads[@]}"
prepare_profile hybrid_mpi_4 "${hybrid_mpi_4_threads[@]}"
prepare_profile hybrid_mpi_16 "${hybrid_mpi_16_threads[@]}"
