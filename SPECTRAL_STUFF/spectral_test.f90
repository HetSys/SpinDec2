!gfortran -std=f2018 -fbounds-check src/Input_parser.f90 src/checkpoint.f90 src/free_energy.f90 src/grid_mod.f90 src/potentials.f90 src/io.f90 src/cahn_hilliard.f90 src/main.f90 `nf-config --fflags --flibs`

program main
    use spectral
    use iso_fortran_env

    implicit none

    real(kind=real64) ,dimension(:,:), allocatable:: c_in,c_out_exp
    real(kind=real64) ,dimension(:,:), allocatable:: c_in_p
    real(kind=real64) ,dimension(:,:), allocatable:: c_out
    real(kind=real64) ,dimension(:), allocatable:: coeffs
    real(kind=real64) :: M,k,dt
    allocate(c_in(3,3))
    allocate(c_in_p(3,3))
    allocate(c_out(3,3))
    allocate(c_out_exp(3,3))
    allocate(coeffs(4))
    coeffs = 2
    c_out_exp = 0
    c_in = 2
    c_in(2,2) = 4
    c_in_p = 2
    c_in_p(2,2) = 5

    M = 0.1
    k = 0.1
    dt = 0.1


    call spectral_method_iter(c_in,c_in_p,coeffs,dt,M,k,c_out,1)

    c_out_exp = reshape((/2.3463434502703331, 45.809290092588583,&
     -7.0049297070283778, 45.809290092588583, -153.91998785683825,&
    45.809290092588583, -7.0049297070283831, 45.809290092588597, &
     2.3463434502703282/),shape(c_out_exp))

    if(maxval(c_out-c_out_exp) > 1e-5) then
        print*, "Test 1 failed"
    else
        print*, "Test 1 passed"
    end if

    call spectral_method_iter(c_in,c_in_p,coeffs,dt,M,k,c_out,0)
    c_out_exp = reshape((/2.4698426892921308, 2.2611576236563615, 2.3736714532770624,&
            2.2611576236563615, 0.93500788690283210, 2.2611576236563611, 2.3736714532770629,&
            2.2611576236563615, 2.4698426892921304/),shape(c_out_exp))

    if(maxval(c_out-c_out_exp) > 1e-5) then
        print*, "Test 2 failed"
    else
        print*, "Test 2 passed"
    end if

end program main
