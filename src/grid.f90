module grid

    use iso_fortran_env

    implicit none

    ! Define module constants
    real(real64), parameter :: pi = 3.1415926535897932

    ! Variables needed for MPI
    real(real64), dimension(:,:), allocatable :: global_grid_conc
    real(real64), dimension(:,:), allocatable :: local_grid_conc
    real(real64), dimension(:,:), allocatable :: conc_halo, Q_halo, M_halo
    
    integer :: grid_domain_size
    integer, dimension(2) :: grid_domain_start, grid_domain_end

contains

    subroutine get_seed(seed_in)
        ! Subroutine that determines random seed used
        ! If seed provided is -1 a random seed is generated
        ! Otherwise the provided seed is used

        implicit none

        integer :: seed_in
        real(real64) :: u
        integer :: n
        integer, dimension(:), allocatable :: seed

        if (seed_in .eq. -1) then
            call random_number(u)
            seed_in = int((999999 - 100000) * u + 100000)
        end if

        ! Set seed
        call random_seed(size=n)
        allocate (seed(n))
        seed = seed_in
        call random_seed(put=seed)

    end subroutine get_seed
    
    subroutine rand_uniform(x,c_min,c_max)

	! Subroutine to get a standard normal random number

	    implicit none

	    real(real64), intent(in) :: c_min, c_max
	    real(real64), intent(out) :: x
	    real(real64) :: u1

	    call random_number(u1)

	    u1 = 1.0 - u1

	    x = (c_max-c_min)*u1 + c_min
 
    end subroutine rand_uniform

    subroutine rand_stdnormal(x)
        ! Subroutine to get a standard normal random number

        implicit none

        real(real64), intent(out) :: x
        real(real64) :: u1, u2

        call random_number(u1)
        call random_number(u2)

        u1 = 1 - u1
        u2 = 1 - u2

        x = sqrt(-2 * log(u1)) * cos(2 * pi * u2)

    end subroutine rand_stdnormal

    subroutine rand_normal(x, mean, std)
        ! Subroutine to get a random number from a normal distribution
        ! given the mean and standard deviation

        implicit none

        real(real64), intent(in) :: mean, std
        real(real64), intent(out) :: x

        call rand_stdnormal(x)
        x = x * std + mean

    end subroutine

    subroutine grid_init(grid, Nx, Ny, c_min, c_max)
        ! Subroutine to initialise the concentration grid
        ! The first input is an array with dimensions Nx x Ny to store concentrations
        ! Nx and Ny are the grid dimensions
        ! c_min is the lower bound for the concentration
        ! c_upper is the upper bound for the concentration
        ! seed_in is a seed for the random number generator

        implicit none

        integer, intent(in) :: Nx, Ny
        real(real64), dimension(Nx, Ny), intent(out) :: grid
        real(real64), intent(in) :: c_min, c_max
        real(real64) :: x ! dummy variable for concnetration
        integer :: i, j ! counters

        do i = 1, Nx
            do j = 1, Ny
                call rand_uniform(x,c_min,c_max)
                grid(i, j) = x
            end do
        end do

    end subroutine grid_init

    ! MPI grid subroutines

    subroutine grid_initialise_local(Nx,Ny,c_min,c_max,p,my_rank,my_rank_coords)

        integer, intent(in) :: Nx, Ny, p
        integer, intent(in) :: my_rank
        integer, dimension(2), intent(in) :: my_rank_coords
        real(real64), intent(in) :: c_min, c_max
        integer :: proot
        integer :: ierr

        proot = int(real(sqrt(real(p,kind=real64)),kind=real64)+0.5)

        grid_domain_size = (Nx*Ny)/proot
        grid_domain_start(:) = 1 + (my_rank_coords(:))*grid_domain_size
        grid_domain_end(:) = grid_domain_start(:) + grid_domain_size - 1

        if (my_rank == 0) then
            print*, "Size of each local processor grid is: ", grid_domain_size, " x ", grid_domain_size
        end if

        ! Set up four subdomains on current rank

        ! Allocate concnentration grid on current rank
        allocate(local_grid_conc(grid_domain_size,grid_domain_size),stat=ierr)
        if (ierr /= 0) then
            print*, "Error: allocating local_grid_con failed on rank ", my_rank
            stop
        end if

        ! Allocate four halo arrays for concentration
        allocate(conc_halo(grid_domain_size,4),stat=ierr)
        if (ierr /= 0) then
            print*, "Error: allocating conc_conc failed on rank ", my_rank
            stop
        end if

        ! Allocate four halo arrays for Q
        allocate(Q_halo(grid_domain_size,4),stat=ierr)
        if (ierr /= 0) then
            print*, "Error: allocating Q_conc failed on rank ", my_rank
            stop
        end if

        ! Allocate four halo arrays for 
        allocate(M_halo(grid_domain_size,4),stat=ierr)
        if (ierr /= 0) then
            print*, "Error: allocating M_halo failed on rank ", my_rank
            stop
        end if

        ! Initiate grid using a uniform distribution on current rank
        call grid_init(global_grid_conc,Nx,Ny,c_min,c_max)

    end subroutine grid_initialise_local

    
    subroutine grid_initialise_global(Nx,Ny,my_rank)

        integer, intent(in) :: Nx, Ny, my_rank
        integer :: ierr

        ! allocate grid array on rank zero
        if (my_rank == 0) then
            allocate(global_grid_conc(Nx,Ny),stat=ierr)
            if(ierr/=0) stop "Error: allocating global_grid_conc failed"
        end if
    end subroutine grid_initialise_global



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
        end if

    end subroutine global_grid_deallocate


    subroutine local_grid_deallocate(my_rank)
        
        integer, intent(in) :: my_rank
        integer :: ierr

        ! Deallocate memory allocated for each local concentration grid
        deallocate(local_grid_conc, stat=ierr)

        if (ierr /= 0) then
            print*, "Error: deallocating local_grid_conc failed on rank ", my_rank
            stop
        end if

    end subroutine local_grid_deallocate

end module grid
