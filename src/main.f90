! Code adapted from PX425 - Assignment 4 (2021) By Dr. David Quigley

program main

    use iso_fortran_env
    use grid_mpi
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
        print*, "A square grid is required. Nx mus equal Ny"
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

    ! Deallocate local and global grids
    call local_grid_deallocate(my_rank)
    call global_grid_deallocate(my_rank)

    ! Finalise MPI
    call comms_finalise()

end program main