! Code adapted from PX425 - Assignment 4 (2021) By Dr. David Quigley

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
    use comms

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
    real(real64) :: stab !Stabilization_Term
    real(real64) :: bfe, df_tol!Placeholder (These were in the input file but df_tol hasn't been used in any code)
    integer :: Nx, Ny, Nt, Nc
    integer :: k, count,thread ! counters
    integer :: cint, random_seed, err, use_input, current_iter, ncerr !checkpointing_interval, random seed,error var
    character(len=128) :: cpi, cpo ! checkpointing files
    character(len=128) :: problem
    integer :: proot ! sqrt of number of processors

    thread =  fftw_init_threads()
    ! Only run files in test for now
    call read_params("input.txt",problem, c_min, c_max, a, nx, &
                     ny, ma, mb,ea,eb,T_min,T_max, kappa, bfe, cint, cpi, cpo, t_end, dt, df_tol,stab, random_seed, use_input, err)

    if (err == -1) then
        print *, "There was an issue with the input file please check and try again"
        stop
    end if

    current_iter = 2
    ! come back to this
    Nt = floor(t_end / dt)
    c0 = (c_min+c_max)/2
    c_std = 0

    ! Allocate T grid
    if (.not. allocated(T)) then
        allocate (T(Nx, Ny))
        T = 0.0
    end if

    if (problem == 'Temp') then
        call grid_init(T, Nx, Ny, T_min, T_max)
    end if

    if (cpi /= "") then
        call read_checkpoint_in(c, mu, T,F_tot, cpi,problem, c0, a, nx, &
                                ny, ma, mb, kappa, bfe, cint, cpo, t_end, dt, df_tol, current_iter, random_seed, use_input, ncerr)
        if (ncerr /= nf90_noerr) then
            print *, "There was an error reading the checkpoint file."
            stop
        end if
    end if

    ! Set seed
    call get_seed(random_seed)

    dx = 1.0_real64/real(nx)
    dy = 1.0_real64/real(ny)

    if( problem /= 'Spectral') then
        !if (dt > min(0.1*dx**4, 0.1*dy**4)) then
        !    print*, 'Warning time-step unstable, setting to default stable value'
        !    dt = min(0.1*dx**4, 0.1*dy**4)
        !end if
    end if

    !Find polynomial coefficients size
    Nc = size(a)
    ! Allocate c grid
    if (.not. allocated(c)) then
        allocate (c(Nx, Ny, Nt))
        c = 0.0
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

    if (problem == 'Temp') then
        call grid_init(T, Nx, Ny, T_min, T_max)
    end if

    ! Check Nx = Ny ie. we have a square grid

    if (Nx /= Ny) then
        print*, "A square grid is required for MPI. Nx must equal Ny"
        stop
    end if

    ! Initialise MPI
    ! Get my_rank and p
    call comms_initialise()

    ! Square root of number of processors
    proot = int(real(sqrt(real(p,kind=real64)),kind=real64)+0.5)

    ! Checks on p and input grid dimensions
    ! Check number of processors is square
    if (proot*proot /= p) then
        if (my_rank == 0) print*, "Error: Number of processors must be exact square"
        stop
    end if

    ! Check that the grid size is divisible by 2*sqrt(p)
    if (mod(Nx*Ny,2*proot) /= 0) then
        if (my_rank == 0) print*, "Error: Number of grid points (Nx*Ny) show be divisible by 2*sqrt(p)"
        stop
    end if

    ! Set up Cartesian communicator
    call comms_processor_map()

    ! write(*,'("MPI rank ",I3," has domain coordinates : ",2I3)')my_rank,my_rank_coords
    ! write(*,'("MPI rank ",I3," has domain neighbours  : ",4I4)')my_rank,my_rank_neighbours

    ! Set up local grids
    call grid_initialise_local(Nx,Ny,p,my_rank,my_rank_coords)

    ! Set up global grid
    ! Initialise global concentration grid
    ! c ~ U(c_min,c_max)
    call grid_initialise_global(Nx,Ny,c_min,c_max,my_rank)

    ! Start calculations

    !Get initial bulk free energy for current rank using local grid
    call bulk_free_energy(f_b,local_grid_conc,a)

    ! Calculate Initial F(t)
    call total_free_energy(local_F, c(:, :, 1), f_b, dx, dy, kappa)

    !Gather total free energy
    call comms_get_global_F(F_tot(1),local_F)

    


    deallocate (f_b)

    count = 0

    ! Grid evolution
    do k = current_iter, Nt

        ! Get bulk chemical potentials
        call bulk_potential(mu,local_grid_conc, a)

        !Store concentration from neighbor ranks in halo_swaps
        call comms_halo_swaps(local_grid_conc,conc_halo)

        ! Get total chemical potentials
        call total_potential(Q, mu,local_grid_conc, dx, dy, Kappa,conc_halo)

        ! Get Mobility Field
        call Mobility(M,MA,MB, EA, EB, c0, local_grid_conc, T, problem)

        !Store Q and M from neighbor ranks 
        call comms_halo_swaps(Q,Q_halo)

        call comms_halo_swaps(M,M_halo)

        ! Get new concentrations for current timesteps
        call time_evolution_new(local_grid_conc,c_new,M,Q,dx,dy,dt, Nx, Ny,Q_halo,M_halo)

        ! set grid to c_new
        local_grid_conc(:, :) = c_new(:, :)

        !Get global concentration grid
        call comms_get_global_grid()

        if (my_rank == 0) then
            c(:,:,k) = global_grid_conc(:,:)
        end if
        
        ! Get Bulk Free Energy over space
        call bulk_free_energy(f_b, c_new, a)

        ! Calculate F(t)
        call total_free_energy(local_F, c_new, f_b, dx, dy, kappa,conc_halo)

        call comms_get_global_F(F_tot(k),local_F)

        deallocate (f_b)

        if (my_rank == 0) then
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
        end if

    end do

    ! Deallocate local and global grids
    call local_grid_deallocate(my_rank)
    call global_grid_deallocate(my_rank)

    ! Finalise MPI
    call comms_finalise()

end program main
