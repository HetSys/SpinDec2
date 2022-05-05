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

        real(real64) :: u
        integer :: seed_in, seed_out

        if (seed_in .EQ. -1) then
        call random_number(u)
            seed_in = int((999999-100000)*u + 100000)
            seed_out = seed_in
        end if

    end subroutine get_seed

    subroutine rand_stdnormal(x,seed_in)
        
        ! Subroutine to get a standard normal random number

        implicit none

        real(real64) :: x
        integer :: seed_in
        real(real64) :: u1, u2
        integer :: n
        integer , dimension(:), allocatable :: seed

        call random_seed(size=n)
        allocate(seed(n))
        seed = seed_in
        call random_seed(put=seed)
        
        call random_number(u1)
        call random_number(u2)

        u1 = 1 - u1
        u2 = 1 - u2
         
        x = sqrt(-2*log(u1))*cos(2*pi*u2)

    end subroutine rand_stdnormal

    subroutine rand_normal(x,seed_in,mean,std)

        ! Subroutine to get a random number from a normal distribution
        ! given the mean and standard deviation

        implicit none

        integer :: seed_in
        real(real64):: x, mean, std

        call rand_stdnormal(x,seed_in)        
        x = x*std + mean

    end subroutine

    subroutine grid_init(grid,Nx,Ny,C,C_std,seed_in)
        ! Subroutine to initialise the concentration grid
        ! The first input is an array with dimensions Nx x Ny to store concentrations
        ! Nx and Ny are the grid dimensions
        ! C is the mean concentration
        ! C_std is the standard deviation of the concentration distribution
        ! seed_in is a seed for the random number generator

        implicit none

        integer, intent(in) :: Nx, Ny
        integer :: seed_in
        real(real64), intent(in) :: C, C_std
        real(real64), dimension(Nx,Ny), intent(out) :: grid
        real(real64) :: x ! dummy variable for concnetration
        integer :: i,j ! counters

        do i = 1, Nx
            do j = 1, Ny
                call rand_normal(x, seed_in, C, C_std)
                grid(i,j) = x
                seed_in = seed_in + 100000
            end do
        end do
        
    end subroutine grid_init
    
end module grid
