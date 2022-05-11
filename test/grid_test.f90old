module grid_test

use iso_fortran_env
use grid

implicit none

contains

subroutine test_rand_seed(seed_in,test_grid_1,test_grid_2,Nx,Ny,C,C_std)

    ! Subroutine to test setting the random seed
    ! Check that setting the seed gives the same
    ! initial grid every time

    real(real64), intent(in) :: C, C_std
    integer, intent(in) :: Nx, Ny
    integer, intent(in) :: seed_in
    real(real64), dimension(Nx,Ny) :: test_grid_1, test_grid_2
    integer :: i,j ! counters
    
    call get_seed(seed_in)

    call grid_init(test_grid_1,Nx,Ny,C,C_std)

    call grid_init(test_grid_2,Nx,Ny,C,C_std)

    ! Compare 2 grids
    do i = 1, Nx
        do j = 1, Ny
            if (abs(test_grid_1(i,j) - test_grid_2(i,j)) < 1e-05) then
                print*, "seed setting passes test"
                exit ! no need to continue check if condition is not satisfied
            else
                print*, "seed setting does not pass test"
            end if 
        end do
    end do

end subroutine test_rand_seed


subroutine test_stdnormal(mean,std)
    ! subroutine to check that rand_stdnormal
    ! gives a standard normal distribution

    real(real64), intent(out) :: mean, std 
    integer :: N = 1000 ! sample size
    real(real64) :: x_i, sum, diff
    integer :: i ! counter

    sum = 0.0
    diff = 0.0

    ! ! Optional - write values to file
    open (99, file = "stdnormal_test.txt")
    do i = 1, N
        call rand_stdnormal(x_i)
        write(99,*) x_i
        sum = sum + x_i
        diff = diff + ((x_i-0.0)**2)
    end do
    close(99)

    ! mean
    mean = sum/N
    ! std
    std = sqrt((diff/N))

    if ((abs(mean-0.0) <= 1e-1) .and. (abs(std-1.0) <= 1e-1)) then
        print*, "rand_stdnormal passes test"
    else
        print*, "rand_stdnormal doesn not pass test"
        ! print*, "mean ", mean
        ! print*, "std ", std
    end if

end subroutine test_stdnormal


subroutine test_rand_normal(mean,std,mean_0,std_0)
    ! subroutine to check that rand_normal
    ! gives a standard normal distribution
    ! given a mean and standard deviation

    real(real64), intent(in) :: mean_0, std_0
    real(real64), intent(out) :: mean, std 
    integer :: N = 1000 ! sample size
    real(real64) :: x_i, sum, diff
    integer :: i ! counter

    sum = 0.0
    diff = 0.0

    ! Optional - write values to file
    open (99, file = "normal_test.txt")
    do i = 1, N
        call rand_normal(x_i,mean_0,std_0)
        write(99,*) x_i
        sum = sum + x_i
        diff = diff + ((x_i-mean_0)**2)
    end do
    close(99)

    ! mean
    mean = sum/N
    ! std
    std = sqrt(diff/N)

    if ((abs(mean-mean_0) <= 1e-3) .and. (abs(std-std_0) <= 1e-3)) then
        print*, "rand_normal passes test"
    else
        print*, "rand_normal doesn not pass test"
        ! print*, "mean ", mean
        ! print*, "std ", std
    end if

end subroutine test_rand_normal


subroutine test_grid_init(grid,Nx,Ny,C_mean,C_std,mean,std)
    ! Subroutine to test the concentration grid
    ! is initialized with the correct distribution

    integer, intent(in) :: Nx, Ny
    real(real64), dimension (Nx,Ny) :: grid
    real(real64), intent(in) :: C_mean, C_std
    real(real64), intent(out) :: mean, std
    real(real64) :: sum, diff
    integer :: i,j ! counters

    ! initialize sum and diff
    sum = 0.0
    diff = 0.0

    ! initialize grid for test
    call grid_init(grid,Nx,Ny,C_mean,C_std)

    do i = 1, Nx
        do j = 1, Ny
            sum = sum + grid(i,j)
            diff = diff + ((grid(i,j)-C_mean)**2)
        end do
    end do

    ! calculate mean grid concentration
    mean = sum/(Nx*Ny)
    ! calculate standard deviation
    std = sqrt(diff/(Nx*Ny))

    if ((abs(mean-C_mean) < 1e-2) .and. (abs(std-C_std) < 1e-2)) then
        print*, "grid_init passes test" 
    else
        print*, "grid_init does not pass test"    
    end if

end subroutine test_grid_init

end module grid_test

program main

    use iso_fortran_env
    use grid
    use grid_test

    implicit none

    real(real64), dimension(:,:), allocatable :: c_grid ! conc. grid
    ! real(real64), dimension(:,:), allocatable :: c_new ! new conc. grid
    ! real(real64), dimension(:,:), allocatable :: mu   ! bulk chem. pot.
    ! real(real64), dimension(:,:), allocatable :: Q    ! total chem. pot.
    ! real(real64), dimension(:,:), allocatable :: dQ   ! 2nd derivative of Q
    real(real64), dimension(:), allocatable :: a
    real(real64) :: C_mean, C_std
    real(real64) :: dx, dy, dt ! spatial and temporal grid spacings
    real(real64) :: Kappa ! free energy gradient parameter
    real(real64) :: M ! Mobility
    integer :: Nx, Ny
    integer :: t_max ! number of timesteps
    integer :: seed_in ! seed for random number generator
    ! integer :: i,j,t ! counters
    ! character(len=100) :: out_file

    ! Testing variables
    real(real64) :: mean, std

    Nx = 10
    Ny = 10

    dx = 0.01
    dy = 0.01
    dt = 1e-12 ! 1 picosecond timestep

    t_max = 10000

    Kappa = 1.6

    M = 1.5

    C_mean = 0.7
    C_std = 0.01

    seed_in = -1

    ! Set seed
    call get_seed(seed_in)

    print*, "Random seed: ", seed_in

    ! Try a 4th order polynomial
    allocate(a(6))
    a = (/5.0,2.1,2.2,2.3,2.4,2.5/)

    ! Allocate grid
    allocate(c_grid(Nx,Ny))

    ! TESTS

    call test_grid_init(c_grid,Nx,Ny,C_mean,C_std,mean,std)

    call test_rand_normal(mean,std,C_mean,C_std)

    call test_stdnormal(mean,std)

end program main