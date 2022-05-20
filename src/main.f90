! Code adapted from PX425 - Assignment 4 (2021) By Dr. David Quigley

program main

    use iso_fortran_env
    use grid
    use comms

    implicit none

    integer :: Nx, Ny
    integer :: seed_in
    real(real64) :: c_min, c_max
    integer :: i,j ! counters
    integer :: proot

    Nx = 4
    Ny = 4

    ! Check Nx = Ny ie. we have a square grid

    if (Nx /= Ny) then
        print*, "A square grid is required. Nx must equal Ny"
        stop
    end if

    seed_in = 123456

    c_min = 0.0
    c_max = 1.0

    call get_seed(seed_in)

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

    write(*,'("MPI rank ",I3," has domain coordinates : ",2I3)')my_rank,my_rank_coords
    write(*,'("MPI rank ",I3," has domain neighbours  : ",4I4)')my_rank,my_rank_neighbours

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
        call time_evoloution_new(local_grid_conc,c_new,M,Q,dx,dy,dt, Nx, Ny,Q_halo,M_halo)

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
