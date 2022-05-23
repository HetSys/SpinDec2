module grid

    use iso_fortran_env

    implicit none

    ! Define module constants
    real(real64), parameter :: pi = 3.1415926535897932

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
        integer, dimension(:),allocatable::seed_o

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

end module grid
