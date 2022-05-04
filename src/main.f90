program main

    use iso_fortran_env
    use cahn_hilliard
    use grid_modules
    use potentials
    use io
    use free_energy

    implicit none

    real(real64), dimension(:,:,:), allocatable :: c ! conc. grid
    real(real64), dimension(:,:), allocatable :: c_new ! new conc. grid
    real(real64), dimension(:,:), allocatable :: mu   ! bulk chem. pot.
    real(real64), dimension(:,:), allocatable :: Q    ! total chem. pot.
    real(real64), dimension(:,:), allocatable :: dQ   ! 2nd derivative of Q
    real(real64), dimension(:), allocatable :: a ! user inputted polynomial coefficients
    real(real64), dimension(:), allocatable :: F_tot ! Total free energy with time
    real(real64), dimension(:,:),allocatable :: f_b
    real(real64) :: c0, c_std !initial grid mean and std sample
    real(real64) :: dx, dy, dt ! spatial and temporal grid spacings
    real(real64) :: Kappa ! free energy gradient parameter
    real(real64) :: t_end !end time
    real(real64) :: M, MA, MB ! Mobility's
    integer :: Nx, Ny, Nt, Nc
    integer :: t_max ! number of timesteps
    integer :: k ! counters


    Nx = 50
    Ny = 50

    dx = 0.01
    dy = 0.01
    dt = 1e-12 ! 1 picosecond timestep
    t_end = 1e-7

    Nt = floor(t_end/dt)

    !t_max = 10000

    Kappa = 1.6

    c0 = 0.7
    c_std = 0.1
    MA = 1
    MB = 1

    !Take constant M to be average estimate of Darken's Equation
    M = (MA*(1-c0) + MB*c0)*c0*(1-c0)



    ! Try a 4th order polynomial
    allocate(a(6))

    a = (/5.0,2.1,2.2,2.3,2.4,2.5/)

    Nc = size(a)

    ! Allocate grid
    allocate(c(Nx,Ny,Nt))

    ! Initialize grid
    call grid_init(c(:,:,1),Nx,Ny,c0,c_std)



    ! Allocate mu
    allocate(mu(Nx,Ny))
    mu = 0.0
    ! Allocate Q
    allocate(Q(Nx,Ny))
    Q = 0.0
    ! Allocate dQ
    allocate(dQ(Nx,Ny))
    dQ = 0.0
    ! Allocate c_new
    allocate(c_new(Nx,Ny))
    c_new = 0.0

    ! Allocate F_tot
    allocate(F_tot(Nt))
    F_tot = 0

    ! Get Initial Bulk Free Energy over space
    call bulk_free_energy(f_b,c(:,:,1),a)

    ! Calculate Initial F(t)
    call total_free_energy(F_tot(1), c(:,:,1), f_b, dx, dy, kappa)

    deallocate(f_b)

    ! Grid evolution
    do k = 2, Nt

        ! Get bulk chemical potentials
        call bulk_potential(mu,c(:,:,k-1),a)

        ! Get total chemical potentials
        call total_potential(Q,mu,c(:,:,k-1),dx,dy,Kappa)

        ! Get second derivative of Q
        call del_Q(dQ,Q,M,dx,dy)

        ! Get new concentrations for current timesteps
        call time_evolution(c(:,:,k-1),c_new,dQ,dt,Nx,Ny)

        ! set grid to c_new
        c(:,:,k) = c_new(:,:)


        ! Get Bulk Free Energy over space
        call bulk_free_energy(f_b,c_new,a)

        ! Calculate F(t)
        call total_free_energy(F_tot(k), c_new, f_b, dx, dy, kappa)

        deallocate(f_b)

    end do


    !Writer for using constant M
    call write_netcdf(c, F_tot, a, Nc, Nx, Ny, Nt, dt, c0, MA, MB, kappa)

end program main
