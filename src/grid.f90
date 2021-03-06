!> Grid module contains subroutines to initialise
!! grid related variables and allocate and deallocate
!! memories for different grids used in the software
module grid

    use iso_fortran_env

    implicit none

    public :: get_seed

    public :: left,right,down,up
    public :: local_grid_conc,grid_domain_size,conc_halo, Q_halo, M_halo
    public :: c,global_grid_conc,grid_domain_start,grid_domain_end

    public :: grid_initialise_local,local_grid_deallocate
    public :: grid_initialise_global,global_grid_deallocate
    public :: Q,M,mu,c_new,T,f_b,write_int

    ! Define module constants
    real(real64), parameter :: pi = 3.1415926535897932

    ! Variables needed for MPI
    real(real64), dimension(:,:), allocatable :: global_grid_conc
    real(real64), dimension(:,:,:), allocatable :: c
    real(real64), dimension(:,:), allocatable :: local_grid_conc
    real(real64), dimension(:,:), allocatable :: Q, M, mu, c_new, T, f_b, c_out
    real(real64), dimension(:,:), allocatable :: conc_halo, Q_halo, M_halo
    integer :: write_int
    ! Constants to define directions
    integer, parameter :: left=1, right=2, down=3, up=4
    ! Constants to define coordinates
    !integer, parameter :: x=1, y=2

    integer :: grid_domain_size
    integer, dimension(2) :: grid_domain_start, grid_domain_end

contains

    !> Subroutine to set random seed based on user input.
    !! If -1 is inputed for seed, a random seed is generated,
    !! otherwise the user provided seed is used.
    !! @param seed_in User inputted seed value
    subroutine get_seed(seed_in)

        integer :: seed_in
        integer :: n
        integer, dimension(:), allocatable :: seed,seed_o

        if (seed_in .eq. -1) then
            call random_seed(size=n)
            allocate (seed(n))
            allocate (seed_o(n))
            !call random_seed(put=seed)
            call random_seed(get=seed_o)
            seed = seed_o(1)
            call random_seed(put=seed)
            !print*, seed
            seed_in = seed(1)
        else
            call random_seed(size=n)
            allocate (seed(n))
            seed = seed_in
            call random_seed(put=seed)
        end if
        ! Set seed


    end subroutine get_seed

    !> Subroutine to output a random number from a uniform distribution
    !! between user provided lower and upper bound values.
    !! @param x Random number output
    !! @param min User inputted lower bound for uniform distribution
    !! @param max User inputted upper bound for uniform distribution
    subroutine rand_uniform(x,min,max)

        real(real64), intent(in) :: min, max
        real(real64), intent(out) :: x
        real(real64) :: u1

        call random_number(u1)

        u1 = 1.0 - u1

        x = (max-min)*u1 + min

    end subroutine rand_uniform

    !> Subroutine to initialise a given grid using values sampled from a uniform
    !! distribution with a given lower and upper bound.
    !! @param grid 2D grid to be populated with random values
    !! @param Nx Grid size in x-direction
    !! @param Ny Grid size in y-direction
    !! @param min User inputted lower bound for uniform distribution
    !! @param max User inputted upper bound for uniform distribution
    !! @param burn (Optional) Parameter to burn initial seed values
    subroutine grid_init(grid, Nx, Ny, min, max, burn)

        ! The first input is an array with dimensions Nx x Ny to store concentrations
        ! Nx and Ny are the grid dimensions
        ! c_min is the lower bound for the concentration
        ! c_upper is the upper bound for the concentration
        ! seed_in is a seed for the random number generator

        implicit none

        integer, intent(in) :: Nx, Ny
        integer, intent(in),optional :: burn
        real(real64), dimension(Nx, Ny), intent(out) :: grid
        real(real64), intent(in) :: min, max
        real(real64) :: x ! dummy variable for concentration
        integer :: i, j ! counters
        if(present(burn)) then
            do i = 1, burn
                call rand_uniform(x,min,max)
            end do
        end if
        do i = 1, Nx
            do j = 1, Ny
                call rand_uniform(x,min,max)
                grid(i, j) = x
            end do
        end do

    end subroutine grid_init

    ! MPI grid subroutines

    !> Subroutine to calculate grid subdomain size,
    !! allocate needed memory for local concentration, Q, and M, and mu grids
    !! and allocate memory for corresponding halos.
    !! @param Nx Grid size in x-direction
    !! @param Ny Grid size in y-direction
    !! @param c_min User inputted miniumum for initial concentration distribution
    !! @param c_max User inputted maximum for initial concentration distribution
    !! @param T_min User inputted miniumum for initial temperature distribution
    !! @param T_max User inputted miniumum for initial temperature distribution
    !! @param problem User inputted problem type
    !! @param p  Number of processors requested for parallel run
    !! @param my_rank Current process' MPI rank
    !! @param my_rank_coords X and Y rank coordinates of the current process
    subroutine grid_initialise_local(Nx,Ny,c_min,c_max,T_min,T_max,problem,p,my_rank,my_rank_coords)

        integer, intent(in) :: Nx, Ny, p
        integer, intent(in) :: my_rank
        integer, dimension(2), intent(in) :: my_rank_coords
        real(real64), intent(in) :: c_min, c_max
        real(real64), intent(in) :: T_min, T_max
        character(len=128), intent(in) :: problem
        integer :: proot
        integer :: ierr

        proot = int(real(sqrt(real(p,kind=real64)),kind=real64)+0.5)

        grid_domain_size = Nx/proot
        grid_domain_start(:) = 1 + (my_rank_coords(:))*grid_domain_size
        grid_domain_end(:) = grid_domain_start(:) + grid_domain_size - 1

        if (my_rank == 0) then
            print*, "Size of each local processor grid is: ", grid_domain_size, " x ", grid_domain_size
        end if

        ! Set up four subdomains on current rank

        ! Allocate local concnentration grid on current rank
        if (.not. allocated(local_grid_conc)) then
            allocate(local_grid_conc(grid_domain_size,grid_domain_size),stat=ierr)
            if (ierr /= 0) then
                print*, "Error: allocating local_grid_conc failed on rank ", my_rank
                stop
            end if
            call grid_init(local_grid_conc,grid_domain_size,grid_domain_size,c_min,c_max,Nx*Ny*my_rank)
        end if

        ! Allocate four halo arrays for concentration
        allocate(conc_halo(grid_domain_size,4),stat=ierr)
        if (ierr /= 0) then
            print*, "Error: allocating conc_halo failed on rank ", my_rank
            stop
        end if

        ! Allocate local c_new grid
        allocate(c_new(grid_domain_size,grid_domain_size),stat=ierr)
        if (ierr /= 0) then
            print*, "Error: allocating c_new failed on rank ", my_rank
            stop
        end if

        ! Allocate local Q grid
        allocate(Q(grid_domain_size,grid_domain_size),stat=ierr)
        if (ierr /= 0) then
            print*, "Error: allocating Q failed on rank ", my_rank
            stop
        end if
        Q = 0.0

        ! Allocate four halo arrays for Q
        allocate(Q_halo(grid_domain_size,4),stat=ierr)
        if (ierr /= 0) then
            print*, "Error: allocating Q_halo failed on rank ", my_rank
            stop
        end if
        Q_halo = 0.0

        ! Allocate local M grid
        allocate(M(grid_domain_size,grid_domain_size),stat=ierr)
        if (ierr /= 0) then
            print*, "Error: allocating M failed on rank ", my_rank
            stop
        end if
        M = 0.0

        ! Allocate four halo arrays for M
        allocate(M_halo(grid_domain_size,4),stat=ierr)
        if (ierr /= 0) then
            print*, "Error: allocating M_halo failed on rank ", my_rank
            stop
        end if
        M_halo = 0.0

        ! Allocate local mu grid
        allocate(mu(grid_domain_size,grid_domain_size),stat=ierr)
        if (ierr /= 0) then
            print*, "Error: allocating mu failed on rank ", my_rank
            stop
        end if
        mu = 0.0

        ! Allocate local T grid
        if (.not. allocated(T)) then
            allocate(T(grid_domain_size,grid_domain_size),stat=ierr)
            if (ierr /= 0) then
                print*, "Error: allocating T failed on rank ", my_rank
                stop
            end if
            T = 0.0
        end if

        ! Allocate local f_b grid
        allocate(f_b(grid_domain_size,grid_domain_size),stat=ierr)
        if (ierr /= 0) then
            print*, "Error: allocating f_b failed on rank ", my_rank
            stop
        end if
        f_b = 0.0

        ! Initialise local temperature grid using a uniform distribution
        ! on current rank
        if (problem == 'Temp') then
            call grid_init(T, grid_domain_size, grid_domain_size, T_min, T_max,Nx*Ny*my_rank+Nx*Ny*p)
        end if

    end subroutine grid_initialise_local

    !> Subroutine to allocate memory for global concentration grid on rank 0.
    !! @param Nx Grid size in x-direction
    !! @param Ny Grid size in y-direction
    !! @param Nt Number of timesteps
    !! @param my_rank Current MPI rank
    subroutine grid_initialise_global(Nx,Ny,Nt,my_rank)

        integer, intent(in) :: Nx, Ny, Nt, my_rank
        integer :: ierr
        ! allocate grid array on rank zero
        if (my_rank == 0) then
            allocate(global_grid_conc(Nx,Ny),stat=ierr)
            if(ierr/=0) stop "Error: allocating global_grid_conc failed"
            if (.not. allocated(c)) then
                allocate(c(Nx,Ny,write_int),stat=ierr)
                if(ierr/=0) stop "Error: allocating c failed"
            end if
        end if
    end subroutine grid_initialise_global

    !> Subroutine to deallocate memory of global concentration grids
    !! @param my_rank Current MPI rank
    subroutine global_grid_deallocate(my_rank)

        integer, intent(in) :: my_rank
        integer :: ierr

        ! Deallocate memory allocated for global grid
        ! on rank 0
        if (my_rank == 0) then
            deallocate(global_grid_conc, stat=ierr)

            if (ierr /= 0) then
                print*, "Error: deallocating global_grid_conc failed on rank 0"
                stop
            end if


            deallocate(c, stat=ierr)

            if (ierr /= 0) then
                print*, "Error: deallocating c failed on rank 0"
                stop
            end if
        end if

    end subroutine global_grid_deallocate

    !> Subroutine to dellocate memory of local grids.
    !! @param my_rank Current MPI rank
    subroutine local_grid_deallocate(my_rank)

        integer, intent(in) :: my_rank
        integer :: ierr

        ! Deallocate memory allocated for each local concentration grid
        deallocate(local_grid_conc, stat=ierr)

        if (ierr /= 0) then
            print*, "Error: deallocating local_grid_conc failed on rank ", my_rank
            stop
        end if

        ! Deallocate memory allocated for each local c_new
        deallocate(c_new, stat=ierr)

        if (ierr /= 0) then
            print*, "Error: deallocating local c_new failed on rank ", my_rank
            stop
        end if

        ! Deallocate memory allocated for each local Q grid
        deallocate(Q, stat=ierr)

        if (ierr /= 0) then
            print*, "Error: deallocating local Q grid failed on rank ", my_rank
            stop
        end if

        ! Deallocate memory allocated for each local M grid
        deallocate(M, stat=ierr)

        if (ierr /= 0) then
            print*, "Error: deallocating local M grid failed on rank ", my_rank
            stop
        end if

        ! Deallocate memory allocated for each local mu grid
        deallocate(mu, stat=ierr)

        if (ierr /= 0) then
            print*, "Error: deallocating local mu grid failed on rank ", my_rank
            stop
        end if

        ! Deallocate memory allocated for each local T grid
        deallocate(T, stat=ierr)

        if (ierr /= 0) then
            print*, "Error: deallocating local T grid failed on rank ", my_rank
            stop
        end if

        ! Deallocate memory allocated for each local f_b grid
        deallocate(f_b, stat=ierr)

        if (ierr /= 0) then
            print*, "Error: deallocating local f_b grid failed on rank ", my_rank
            stop
        end if

        ! Deallocate halo memories
        deallocate(conc_halo, stat=ierr)

        if (ierr /= 0) then
            print*, "Error: deallocating conc_halo failed on rank ", my_rank
            stop
        end if

        deallocate(Q_halo, stat=ierr)

        if (ierr /= 0) then
            print*, "Error: deallocating Q_halo failed on rank ", my_rank
            stop
        end if

        deallocate(M_halo, stat=ierr)

        if (ierr /= 0) then
            print*, "Error: deallocating M_halo failed on rank ", my_rank
            stop
        end if

    end subroutine local_grid_deallocate

end module grid
