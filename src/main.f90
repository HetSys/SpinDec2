! Adapted from code originally written by Dr. David Quigley

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
    real(real64), dimension(:, :), allocatable :: mu_check   ! bulk chem. pot.
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
    real(real64) :: bfe!Placeholder (These were in the input file but df_tol hasn't been used in any code)
    integer :: Nx, Ny, Nt, Nc,no_threads,write_num,write_freq_c,last_write
    integer :: k, count,thread ,write_count,write_prev,write_next! counters
    integer :: cint, random_seed, err, use_input, current_iter, ncerr,singl, write_freq !checkpointing_interval, random seed,error var
    character(len=128) :: cpi, cpo ! checkpointing files
    character(len=128) :: problem
    integer :: proot, file_id,i ! sqrt of number of processors
    real(real64) :: t1,t2,w1,w2

    ! Initialise MPI
    ! Get my_rank and p
    !no_threads = omp_get_num_threads()
    call comms_initialise()

    t1 = MPI_Wtime()


    thread =  fftw_init_threads()

    ! Only run files in test for now
    call read_params("input.txt",problem, c_min, c_max, a, nx, &
                     ny, ma, mb,ea,eb,T_min,T_max, kappa, bfe, cint, &
                     cpi, cpo, t_end, dt,stab, random_seed, use_input,singl,write_freq,err)

    if (err == -1) then
        print *, "There was an issue with the input file please check and try again"
        call comms_finalise()
        stop
    end if

    current_iter = 2
    ! come back to this
    write_int = 2

    if(my_rank == 0) then
        write_count = 1
        write_freq_c = write_freq
        last_write=1
    end if

    if(problem == "spectral") then
        !call omp_set_num_threads(no_threads)
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
                                ny, ma, mb, kappa, bfe, cint, cpo, t_end, dt, &
                                current_iter, random_seed, use_input, &
                                last_write,write_freq,write_freq_c,singl&
                                ,write_count,stab, ncerr)
        if (ncerr /= nf90_noerr) then
            print *, "There was an error reading the checkpoint file."
            call comms_finalise()
            stop
        end if
        if(my_rank == 0) then
            write_count = write_count +1
            if(write_count > write_int) then
                write_count = 1
            end if
            write_freq_c = write_freq_c +1
        end if
    end if

    Nt = floor(t_end / dt)
    c0 = (c_min+c_max)/2




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


    if (cpi /= "") then
        call read_checkpoint_data(c, T,F_tot, cpi,use_input,ncerr)
        if (ncerr /= nf90_noerr) then
            print *, "There was an error reading the checkpoint file."
            call comms_finalise()
            stop
        end if
    end if

    call comms_get_global_grid()


    if (my_rank==0) then
        c(:,:,write_count) = global_grid_conc(:,:)
    end if


    if (cpi == "") then
        !Get initial bulk free energy for current rank using local grid
        call bulk_free_energy(f_b,local_grid_conc,a)

        !Store concentration from neighbor ranks in halo_swaps
        call comms_halo_swaps(local_grid_conc,conc_halo)

        ! Calculate Initial F(t)
        call total_free_energy(local_F, local_grid_conc,f_b, dx, dy, kappa,conc_halo)

        !Gather total free energy
        call comms_get_global_F(local_F,global_F)

        if (my_rank == 0) then
            F_tot(1) = global_F
            print*,global_F,F_tot(1:1)
        end if

        if (my_rank == 0) then
            if(write_freq == 0) then
                write_num = 1
            else
                if(mod(Nt,write_freq) == 0) then
                    write_num = Nt/write_freq
                else
                    write_num = Nt/write_freq +1
                end if
            end if
            call write_netcdf_parts_setup(c, F_tot, a, Nc, Nx, Ny, write_num, dt, c0, MA, MB, kappa,t_end,file_id)
            !call write_netcdf_parts(c(:,:,write_count:write_count), F_tot(1:1),1,file_id)
        end if
    else
        if (my_rank == 0) then
            call open_netcdf(file_id)
        end if
    end if
    count = 0
    ! Grid evolution
    do k = current_iter, Nt
        if(write_count == write_int) then
            write_next = 1
        else
            write_next = write_count+1
        end if
        if(write_count == 1) then
            write_prev = write_int
        else
            write_prev = write_count-1
        end if
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
        if(problem /= "spectral") then
            call time_evolution_new(local_grid_conc,c_new,M,Q,dx,dy,dt,grid_domain_size,grid_domain_size,Q_halo,M_halo)

            ! set grid to c_new
            local_grid_conc(:, :) = c_new(:, :)

            call comms_get_global_grid()


            if (my_rank == 0) then
                c(:,:,write_next) = global_grid_conc(:,:)
            end if
        else

            if(k == 2) then
                call spectral_method_iter(c(:, :, write_count),c(:, :,write_count),a,dt,M(1,1),Kappa,c_new,1,stab)
                c(:,:,write_next) = c_new(:,:)
            else

                call spectral_method_iter(c(:, :, write_count),c(:, :,write_prev),a,dt,M(1,1),Kappa,c_new,0,stab)
                c(:,:,write_next) = c_new(:,:)

                !print*, c_new(:,:)
            end if
        end if

        if(p == 1) then
            local_grid_conc(:,:) = c(:,:,write_next)
        end if



        !Get global concentration grid


        ! Get Bulk Free Energy over space
        call bulk_free_energy(f_b,local_grid_conc, a)

        ! Calculate F(t)
        call total_free_energy(local_F, local_grid_conc, f_b, dx, dy, kappa,conc_halo)

        call comms_get_global_F(local_F,global_F)

        if (my_rank == 0) then
            F_tot(k) = global_F
        end if

        if (my_rank == 0) then
            if(write_count == write_int-1 .and. write_freq == 1) then
                w1 = MPI_Wtime()
                if(singl == 1) then
                    c = sngl(c)
                end if
                call write_netcdf_parts(c(:,:,:), F_tot(k-write_int+1:k),k-write_int+1,file_id)
                w2 = MPI_Wtime()

                if (my_rank == 0) then
                    print *, "Write took", w2-w1
                end if
            end if
            if(write_freq > 1 .and. write_freq==write_freq_c) then
                print*, "Writing timestep at iter", k
                write_freq_c = 0
                call write_netcdf_parts(c(:,:,write_next:write_next), F_tot(k:k),last_write,file_id)
                last_write = last_write+1
            end if
        end if

        if (my_rank == 0) then
            if (count >= cint) then
                call write_checkpoint_file(c, T,F_tot,problem, a, cpo, c0, &
                                           nx, ny, ma, mb, kappa, bfe, Cint, t_end, dt, k, &
                                           random_seed, last_write,write_freq,write_freq_c,singl&
                                           ,write_count,stab,ncerr)

                if (ncerr /= nf90_noerr) then
                    print *, "There was an error writing the checkpoint file."
                    call comms_finalise()
                    stop
                end if

                count = 0
            end if
            if(my_rank == 0) then
                write_count = write_count +1
                if(write_count > write_int) then
                    write_count = 1
                end if
                write_freq_c = write_freq_c +1
            end if
            count = count + 1
        end if

    end do
    if (my_rank == 0) then
        if(write_count /= write_int .and. write_freq == 1) then
            if(singl == 1) then
                c = sngl(c)
            end if
            call write_netcdf_parts(c(:,:,:write_count), &
                                F_tot(Nt-write_count:),Nt-write_count,file_id)

        end if
        if(write_freq == 0) then
            call write_netcdf_parts(c(:,:,write_count:write_count), F_tot(Nt:Nt),1,file_id)
        end if
    end if


    ! Deallocate local and global grids
    call local_grid_deallocate(my_rank)
    call global_grid_deallocate(my_rank)

    t2 = MPI_Wtime()

    if (my_rank == 0) then
        !call close_netcdf(file_id)
        print *, "Calculation took", t2-t1, "seconds on", p, "MPI threads"
    end if

    ! Finalise MPI
    call comms_finalise()

end program main
