program main

    use iso_fortran_env
    use cahn_hilliard
    use grid
    use potentials
    use io
    use free_energy
    use input_params
    use checkpointing

    implicit none

    real(real64), dimension(:, :, :), allocatable :: c ! conc. grid
    real(real64), dimension(:, :), allocatable :: c_new ! new conc. grid
    real(real64), dimension(:, :), allocatable :: mu   ! bulk chem. pot.
    real(real64), dimension(:, :), allocatable :: Q    ! total chem. pot.
    real(real64), dimension(:, :), allocatable :: dQ   ! 2nd derivative of Q
    real(real64), dimension(:), allocatable :: a ! user inputted polynomial coefficients
    real(real64), dimension(:), allocatable :: F_tot ! Total free energy with time
    real(real64), dimension(:, :), allocatable :: f_b
    real(real64) :: c0, c_std !initial grid mean and std sample
    real(real64) :: dx, dy, dt ! spatial and temporal grid spacings
    real(real64) :: Kappa ! free energy gradient parameter
    real(real64) :: t_end !end time
    real(real64) :: M, MA, MB ! Mobility's
    real(real64) :: bfe, df_tol!PLaceholder (These were in the input file but df_tol hasn't been used in any code)
    integer :: Nx, Ny, Nt, Nc
    integer :: k, count ! counters
    integer :: cint, random_seed, err, use_input, current_iter, ncerr !checkpointing_interval, random seed,error var
    character(len=128) :: cpi, cpo ! checkpointing files

    ! Only run files in test for now
    call read_params("../test/input_test.txt", c0, c_std, a, nx, &
                     ny, ma, mb, kappa, bfe, cint, cpi, cpo, t_end, dt, df_tol, random_seed, use_input, err)

    if (err == -1) then
        print *, "There was an issue with the input file please check and try again"
        stop
    end if

    current_iter = 2
    Nt = floor(t_end / dt)

    if (cpi /= "") then
        call read_checkpoint_in(c, mu, F_tot, cpi, c0, c_std, a, nx, &
                                ny, ma, mb, kappa, bfe, cint, cpo, t_end, dt, df_tol, current_iter, random_seed, use_input, ncerr)
        if (ncerr /= nf90_noerr) then
            print *, "There was an error reading the checkpoint file."
            stop
        end if
    end if

    ! Set seed
    call get_seed(random_seed)

    dx = 0.01
    dy = 0.01

    !Take constant M to be average estimate of Darken's Equation
    M = (MA * (1 - c0) + MB * c0) * c0 * (1 - c0)

    !Find polynomial coefficients size
    Nc = size(a)

    ! Allocate grid
    if (.not. allocated(c)) then
        allocate (c(Nx, Ny, Nt))
        c = 0
    end if

    ! Allocate mu
    if (.not. allocated(mu)) then
        allocate (mu(Nx, Ny))
        mu = 0.0
    end if

    ! Allocate F_tot
    if (.not. allocated(F_tot)) then
        allocate (F_tot(Nt))
        F_tot = 0
    end if

    ! Allocate Q
    allocate (Q(Nx, Ny))
    Q = 0.0

    ! Allocate dQ
    allocate (dQ(Nx, Ny))
    dQ = 0.0

    ! Allocate c_new
    allocate (c_new(Nx, Ny))
    c_new = 0.0

    ! Initialize grid
    call grid_init(c(:, :, 1), Nx, Ny, c0, c_std)

    ! Get Initial Bulk Free Energy over space
    call bulk_free_energy(f_b, c(:, :, 1), a)

    ! Calculate Initial F(t)
    call total_free_energy(F_tot(1), c(:, :, 1), f_b, dx, dy, kappa)

    deallocate (f_b)

    count = 0

    ! Grid evolution
    do k = current_iter, Nt

        ! Get bulk chemical potentials
        call bulk_potential(mu, c(:, :, k - 1), a)

        ! Get total chemical potentials
        call total_potential(Q, mu, c(:, :, k - 1), dx, dy, Kappa)

        ! Get second derivative of Q
        call del_Q(dQ, Q, M, dx, dy, Nx, Ny)

        ! Get new concentrations for current timesteps
        call time_evolution(c(:, :, k - 1), c_new, dQ, dt, Nx, Ny)

        ! set grid to c_new
        c(:, :, k) = c_new(:, :)

        ! Get Bulk Free Energy over space
        call bulk_free_energy(f_b, c_new, a)

        ! Calculate F(t)
        call total_free_energy(F_tot(k), c_new, f_b, dx, dy, kappa)

        deallocate (f_b)

        if (count >= cint) then
            call write_checkpoint_file(c, mu, F_tot, a, cpo, c0, c_std, &
                                       nx, ny, ma, mb, kappa, bfe, Cint, t_end, dt, k, df_tol, &
                                       random_seed, ncerr)

            if (ncerr /= nf90_noerr) then
                print *, "There was an error writing the checkpoint file."
                stop
            end if

            count = 0
        end if

        count = count + 1

    end do

    !Writer for using constant M
    call write_netcdf(c, F_tot, a, Nc, Nx, Ny, Nt, dt, c0, MA, MB, kappa)

end program main
