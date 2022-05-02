module grid_modules

    use iso_fortran_env

    implicit none

    ! Define module constants
    real(real64), parameter :: pi = 3.1415926535897932

    contains

    subroutine rand_stdnormal(x)
        
        ! Subroutine to get a standard normal random number

        implicit none

        real(real64) :: x
        real(real64) :: u1, u2

        call random_number(u1)
        call random_number(u2)

        u1 = 1 - u1
        u2 = 1 - u2
         
        x = sqrt(-2*log(u1))*cos(2*pi*u2)

    end subroutine rand_stdnormal

    subroutine rand_normal(x,mean,std)

        ! Subroutine to get a random number from a normal distribution
        ! given the mean and standard deviation

        implicit none

        real(real64):: x, mean, std

        call rand_stdnormal(x)        
        x = x*std + mean

    end subroutine

    subroutine grid_init(grid,Nx,Ny,C,C_std)
        ! Subroutine to initialise the concentration grid
        ! The first input is an array with dimensions Nx x Ny to store concentrations
        ! C is the mean concentration
        ! C_std is the standard deviation of the distribution

        implicit none

        integer :: Nx, Ny
        real(real64), dimension(Nx,Ny) :: grid
        real(real64) :: x ! dummy variable for concnetration
        real(real64) :: C, C_std
        integer :: i,j ! counters

        do i = 1, Nx
            do j = 1, Ny
                call rand_normal(x, C, C_std)
                grid(i,j) = x    
            end do
        end do
        
    end subroutine grid_init

end module grid_modules
