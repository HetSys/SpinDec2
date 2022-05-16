!gfortran -std=f2018 -fbounds-check src/Input_parser.f90 src/checkpoint.f90 src/free_energy.f90 src/grid_mod.f90 src/potentials.f90 src/io.f90 src/cahn_hilliard.f90 src/main.f90 `nf-config --fflags --flibs`

program main
    use spectral
    use iso_fortran_env

    implicit none

    real(kind=real64) ,dimension(:,:), allocatable:: c_in
    real(kind=real64) ,dimension(:,:), allocatable:: c_in_p
    real(kind=real64) ,dimension(:,:), allocatable:: c_out
    real(kind=real64) ,dimension(:), allocatable:: coeffs
    real(kind=real64) :: M,k,dt
    allocate(c_in(3,3))
    allocate(c_in_p(3,3))
    allocate(c_out(3,3))
    allocate(coeffs(4))
    coeffs = 2

    c_in = 2
    c_in(2,2) = 4
    c_in_p = 2
    c_in_p(2,2) = 5

    M = 0.1
    k = 0.1
    dt = 0.1


    call spectral_method_iter(c_in,c_in_p,coeffs,dt,M,k,c_out,1)

    print*,c_out
end program main
