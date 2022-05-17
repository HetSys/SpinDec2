program main

    use iso_fortran_env
    use cahn_hilliard
    use grid
    use potentials
    use spectral
    use io
    use free_energy
    use input_params
    use checkpointing

    implicit none

    real(real64), dimension(:, :, :), allocatable :: c ! conc. grid
    real(real64), dimension(:, :), allocatable :: c_new, c_out ! new conc. grid
    real(real64), dimension(:, :), allocatable :: mu   ! bulk chem. pot.
    real(real64), dimension(:, :), allocatable :: Q    ! total chem. pot.
    real(real64), dimension(:, :), allocatable :: dQ   ! 2nd derivative of Q
    real(real64), dimension(:, :), allocatable :: M    ! Mobility field
    real(real64), dimension(:, :), allocatable :: T ! Temp
    real(real64), dimension(:), allocatable :: a ! user inputted polynomial coefficients
    real(real64), dimension(:), allocatable :: F_tot ! Total free energy with time
    real(real64), dimension(:, :), allocatable :: f_b
    real(real64) :: c0, c_std !initial grid mean and std sample - if using normal dist
    real(real64) :: c_min, c_max !initial grid lower and upper boubds - for default uniform dist
    real(real64) :: T_min, T_max !initial temp grid mean and std sample
    real(real64) :: dx, dy, dt ! spatial and temporal grid spacings
    real(real64) :: Kappa ! free energy gradient parameter
    real(real64) :: t_end !end time
    real(real64) :: MA, MB ! Mobility's
    real(real64) :: EA, EB ! exciation energy
    real(real64) :: bfe, df_tol!Placeholder (These were in the input file but df_tol hasn't been used in any code)
    integer :: Nx, Ny, Nt, Nc
    integer :: k, count,thread ! counters
    integer :: cint, random_seed, err, use_input, current_iter, ncerr !checkpointing_interval, random seed,error var
    character(len=128) :: cpi, cpo ! checkpointing files
    character(len=*), parameter :: problem = "Constant"


    thread =  fftw_init_threads()
    ! Only run files in test for now
    call read_params("input.txt", c0, c_std, a, nx, &
                     ny, ma, mb, kappa, bfe, cint, cpi, cpo, t_end, dt, df_tol, random_seed, use_input, err)

    if (err == -1) then
        print *, "There was an issue with the input file please check and try again"
        stop
    end if

    current_iter = 2
    ! come back to this
    ! Nt = floor(t_end / dt)

    Nt = 1e4

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

    if( problem /= 'Spectral') then
        if (dt > min(0.1*dx**4, 0.1*dy**4)) then
            print*, 'Warning time-step unstable, setting to default stable value'
            dt = min(0.1*dx**4, 0.1*dy**4)
        end if
    end if

    EA = 1.0
    EB = 1.0

    c_min = 0.1
    c_max = 0.9

    T_min = 945
    T_max = 955


    !Find polynomial coefficients size
    Nc = size(a)
    ! Allocate c grid
    if (.not. allocated(c)) then
        allocate (c(Nx, Ny, Nt))
        c = 0.0
    end if

    ! Allocate T grid
    if (.not. allocated(T)) then
        allocate (T(Nx, Ny))
        T = 0.0
    end if

     ! Allocate M grid
    if (.not. allocated(M)) then
        allocate (M(Nx, Ny))
        M = 0.0
    end if

    ! Allocate mu grid
    if (.not. allocated(mu)) then
        allocate (mu(Nx, Ny))
        mu = 0.0
    end if

    ! Allocate F_tot
    if (.not. allocated(F_tot)) then
        allocate (F_tot(Nt))
        F_tot = 0.0
    end if

     ! Allocate Q grid
    if (.not. allocated(Q)) then
        allocate (Q(Nx, Ny))
        Q = 0.0
    end if

     ! Allocate c_new grid
    if (.not. allocated(c_new)) then
        allocate (c_new(Nx, Ny))
        c_new = 0.0
    end if

    if (.not. allocated(c_out)) then
        allocate (c_out(Nx, Ny))
        c_out = 0.0
    end if

    ! Initialize c grid
    call grid_init(c(:, :, 1), Nx, Ny, c_min, c_max)

    if (problem == 'Temp') then
        call grid_init(T, Nx, Ny, T_min, T_max)
    end if

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

        ! Get Mobility Field
        call Mobility(M,MA,MB, EA, EB, c0, c(:, :, k-1), T, problem)

        !print*, M(1,1), M(1,2), M(6,7)

        ! Get new concentrations for current timesteps

        if (problem == 'Spectral') then
            if(k == 2) then
                call spectral_method_iter(c(:, :, 1),c(:, :, 1),a,dt,M(1,1),Kappa,c_out,1)
                c(:,:,2) = c_out
            else
                call spectral_method_iter(c(:, :, k-1),c(:, :, k-2),a,dt,M(1,1),Kappa,c_out,0)
                c(:,:,k) = c_out
            end if
        else
            call time_evoloution_new(c(:, :, k-1),c_new,M,Q,dx,dy,dt, Nx, Ny)

            c(:, :, k) = c_new(:, :)
        end if

        ! Get Bulk Free Energy over space
        call bulk_free_energy(f_b, c(:, :, k), a)

        ! Calculate F(t)
        call total_free_energy(F_tot(k), c(:, :, k), f_b, dx, dy, kappa)

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
