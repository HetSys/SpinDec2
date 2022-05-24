! Code adapted from PX425 - Assignment 4 (2021) By Dr. David Quigley

program main

    use iso_fortran_env
    use cahn_hilliard
    use checkpointing
    use comms
    use free_energy
    use grid
    use input_params
    use io
    use potentials
    use mpi
    use omp_lib
    use spectral

    implicit none

    real(real64), dimension(:, :, :), allocatable :: c_check ! conc. grid
    !real(real64), dimension(:, :), allocatable :: c_new, c_out ! new conc. grid
    real(real64), dimension(:, :), allocatable :: mu_check ,phold  ! bulk chem. pot.
    real(real64), dimension(:, :), allocatable :: Q_check    ! total chem. pot.
    real(real64), dimension(:, :), allocatable :: M_check    ! Mobility field
    real(real64), dimension(:, :), allocatable :: T_check ! Temp
    real(real64), dimension(:), allocatable :: a ! user inputted polynomial coefficients
    real(real64), dimension(:), allocatable :: F_tot ! Total free energy with time
    !real(real64), dimension(:, :), allocatable :: f_b
    real(real64) :: c0, c_std !initial grid mean and std sample - if using normal dist
    real(real64) :: c_min, c_max !initial grid lower and upper boubds - for default uniform dist
    real(real64) :: T_min, T_max !initial temp grid mean and std sample
    real(real64) :: dx, dy, dt ! spatial and temporal grid spacings
    real(real64) :: Kappa ! free energy gradient parameter
    real(real64) :: t_end !end time
    real(real64) :: MA, MB ! Mobility's
    real(real64) :: EA, EB ! exciation energy
    real(real64) :: local_F,global_F
    real(real64) :: stab !Stabilization_Term
    real(real64) :: bfe, df_tol!Placeholder (These were in the input file but df_tol hasn't been used in any code)
    integer :: Nx, Ny, Nt, Nc,no_threads
    integer :: k, count,thread ! counters
    integer :: cint, random_seed, err, use_input, current_iter, ncerr !checkpointing_interval, random seed,error var
    character(len=128) :: cpi, cpo ! checkpointing files
    character(len=128) :: problem
    integer :: proot, file_id ! sqrt of number of processors
    real(real64) :: t1,t2

    ! Initialise MPI
    ! Get my_rank and p
    !no_threads = omp_get_num_threads()
    call comms_initialise()

    t1 = MPI_Wtime()


    thread =  fftw_init_threads()
    ! Only run files in test for now
    call read_params("input.txt",problem, c_min, c_max, a, nx, &
                     ny, ma, mb,ea,eb,T_min,T_max, kappa, bfe, cint, cpi, cpo, t_end, dt, df_tol,stab, random_seed, use_input, err)

    if (err == -1) then
        print *, "There was an issue with the input file please check and try again"
        call comms_finalise()
        stop
    end if

    current_iter = 2
    ! come back to this
    Nt = floor(t_end / dt)
    c0 = (c_min+c_max)/2


    if(problem == "Spectral") then
        call omp_set_num_threads(no_threads)
        if(p > 1) then
            if(my_rank == 0) then
                print*, "Spectral can oly be run with pure OpenMP"
            end if
            call comms_finalise()
            stop
        end if
    end if

    if (cpi /= "") then
        call read_checkpoint_metadata(cpi,problem, c0, a, nx, &
                                ny, ma, mb, kappa, bfe, cint, cpo, t_end, dt, df_tol, current_iter, random_seed, use_input, ncerr)
        if (ncerr /= nf90_noerr) then
            print *, "There was an error reading the checkpoint file."
            call comms_finalise()
            stop
        end if
    end if



    ! Set seed
    call get_seed(random_seed)

    dx = 1.0_real64/real(nx)
    dy = 1.0_real64/real(ny)



    !Find polynomial coefficients size
    Nc = size(a)

    ! Allocate F_tot on rank 0
    if (.not. allocated(F_tot)) then
        allocate (F_tot(Nt))
        F_tot = 0.0
    end if



    ! Check Nx = Ny ie. we have a square grid

    if (Nx /= Ny) then
        print*, "A square grid is required for MPI. Nx must equal Ny"
        call comms_finalise()
        stop
    end if


    ! Square root of number of processors
    proot = int(real(sqrt(real(p,kind=real64)),kind=real64)+0.5)

    ! Checks on p and input grid dimensions
    ! Check number of processors is square
    if (proot*proot /= p) then
        if (my_rank == 0) print*, "Error: Number of processors must be exact square"
        call comms_finalise()
        stop
    end if

    ! Check that the grid size is divisible by 2*sqrt(p)
    if (mod(Nx*Ny,2*proot) /= 0) then
        if (my_rank == 0) print*, "Error: Number of grid points (Nx*Ny) show be divisible by 2*sqrt(p)"
        call comms_finalise()
        stop
    end if


    ! Set up Cartesian communicator
    call comms_processor_map()


    ! Set up local grids
    call grid_initialise_local(Nx,Ny,c_min,c_max,T_min,T_max,problem,p,my_rank,my_rank_coords)

    ! Set up global grid
    ! Initialise global concentration grid
    ! c ~ U(c_min,c_max)


    call grid_initialise_global(Nx,Ny,Nt,my_rank)

    ! Start calculations

    if(.not. allocated(phold)) then
        allocate(phold(1,1))
        phold = 1
    end if

    if (cpi /= "") then
        call read_checkpoint_data(c, phold, T,F_tot, cpi,use_input,ncerr)
        if (ncerr /= nf90_noerr) then
            print *, "There was an error reading the checkpoint file."
            call comms_finalise()
            stop
        end if
    end if

    call comms_get_global_grid()


    if (my_rank==0) then
        c(:,:,2) = global_grid_conc(:,:)
    end if

    if (cpi == "") then
        !Get initial bulk free energy for current rank using local grid
        call bulk_free_energy(local_grid_conc,a)

        !Store concentration from neighbor ranks in halo_swaps
        call comms_halo_swaps(local_grid_conc,conc_halo)

        ! Calculate Initial F(t)
        call total_free_energy(local_F, local_grid_conc, dx, dy, kappa,conc_halo)

        !Gather total free energy
        call comms_get_global_F(local_F,global_F)

        if (my_rank == 0) then
            F_tot(1) = global_F
        end if
        if (my_rank == 0) then
            call write_netcdf_parts_setup(c, F_tot, a, Nc, Nx, Ny, Nt, dt, c0, MA, MB, kappa,file_id)
            call write_netcdf_parts(c(:,:,2), F_tot(1),1,file_id)
        end if
    end if
    count = 0
    ! Grid evolution
    do k = current_iter, Nt
        !print*, k
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
        if(problem /= "Spectral") then
            call time_evolution_new(local_grid_conc,c_new,M,Q,dx,dy,dt,grid_domain_size,grid_domain_size,Q_halo,M_halo)

            ! set grid to c_new
            local_grid_conc(:, :) = c_new(:, :)

            call comms_get_global_grid()


            if (my_rank == 0) then
                c(:,:,3) = global_grid_conc(:,:)
            end if
        else
            if(k == 2) then
                call spectral_method_iter(c(:, :, 2),c(:, :, 2),a,dt,M(1,1),Kappa,c_new,1,stab)
                c(:,:,3) = c_new(:,:)
            else
                call spectral_method_iter(c(:, :, 2),c(:, :, 1),a,dt,M(1,1),Kappa,c_new,0,stab)
                c(:,:,3) = c_new(:,:)
            end if
        end if

        if(p == 1) then
            local_grid_conc(:,:) = c(:,:,3)
        end if



        !Get global concentration grid


        ! Get Bulk Free Energy over space
        call bulk_free_energy(local_grid_conc, a)

        ! Calculate F(t)
        call total_free_energy(local_F, local_grid_conc, dx, dy, kappa,conc_halo)

        call comms_get_global_F(local_F,global_F)

        if (my_rank == 0) then
            F_tot(k) = global_F
        end if

        if (my_rank == 0) then
            call write_netcdf_parts(c(:,:,3), F_tot(k),k,file_id)

            c(:,:,1) = c(:,:,2)
            c(:,:,2) = c(:,:,3)
        end if

        if (my_rank == 0) then
            if (count >= cint) then
                call write_checkpoint_file(c, phold, T,F_tot,problem, a, cpo, c0, &
                                           nx, ny, ma, mb, kappa, bfe, Cint, t_end, dt, k, df_tol, &
                                           random_seed, ncerr)

                if (ncerr /= nf90_noerr) then
                    print *, "There was an error writing the checkpoint file."
                    call comms_finalise()
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
    deallocate(phold)

    t2 = MPI_Wtime()

    if (my_rank == 0) then
        call close_netcdf(file_id)
        print *, "Calculation took", t2-t1, "seconds on", p, "MPI threads"
    end if

    ! Finalise MPI
    call comms_finalise()

end program main
