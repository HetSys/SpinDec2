program main
    use input_params
    use checkpointing
    use netcdf
    use iso_fortran_env
    use input_tester
    use checkpoint_tester
    use test_potentials
    use test_free_energy
    use test_derivatives
    use grid_test
    use ch_test

    implicit none

    real(real64), dimension(3) :: a_1
    real(real64), dimension(3, 3) :: c_1
    real(real64), dimension(3, 3) :: expected_bulk_1, expected_total_1, expected_fb_1
    real(real64) :: dx, dy, kappa, dt
    real(real64), dimension(5) :: a_2
    real(real64), dimension(4, 4) :: c_2
    real(real64), dimension(4, 4) :: expected_bulk_2, expected_total_2, expected_fb_2
    real(real64), dimension(3) :: a_3
    real(real64), dimension(3,3)  :: c_3,expected_bulk_3,expected_total_3
    real(real64) :: expected_total_F_1, expected_total_F_2
    integer :: test_num
    real(real64), dimension(:,:), allocatable :: c_new ! new conc. grid
    real(real64), dimension(:,:), allocatable :: dQ   ! 2nd derivative of Q
    real(real64), dimension(:), allocatable :: a
    real(real64) :: C_mean, C_std
    real(real64) :: M ! Mobility
    integer :: Nx, Ny
    integer :: t_max ! number of timesteps
    real(real64), dimension(:,:), allocatable :: c_grid ! conc. grid
    integer :: seed_in ! seed for random number generator
    ! Testing variables
    real(real64) :: mean, std
    real(real64), dimension(:,:), allocatable :: Q_test, dQ_expected
    real(real64), dimension(:,:), allocatable :: c_test, c_expected
    integer :: i,j ! counters
    print*,'########################################################'
    print *, "Starting input testing"
    print*,'--------------------------------------------------------'
    call input_test()

    print*,'########################################################'
    print *, "Starting checkpoint testing"
    print*,'--------------------------------------------------------'
    call checkpoint_test()

    print*,'########################################################'
    print *, 'Starting potentials unit testing.'
    print*,'--------------------------------------------------------'

    !Test 1 for bulk potential
    a_1 = (/1.0, 2.0, 3.0/)
    c_1 = reshape((/1.0, 2.0, 3.0, 2.0, 3.0, 6.0, 1.0, 2.0, 4.0/), shape(c_1))
    expected_bulk_1 = reshape((/8.0, 14.0, 20.0, 14.0, 20.0, 38.0, 8.0, 14.0, 26.0/), shape(c_1))
    test_num = 1
    call test_bulk_potential(c_1, a_1, expected_bulk_1, test_num)

    !Test 2 for bulk potential
    a_2 = (/0.2, 0.5, 0.3, 1.4, 1.6/)
    c_2 = reshape((/0.16, 0.23, 0.32, 0.45, 0.35, 0.42, 0.71, 0.12, 0.24, 0.45, 0.64, 0.83, 0.22, 0.42, 0.18, 0.71/), shape(c_2))
    expected_bulk_2 = reshape((/0.7297, 0.9380, 1.3317, 2.2037, 1.4989, 1.9670, 5.3338, 0.6435, 0.9743, 2.2037, &
                                4.2820, 7.5508, 0.9034, 1.9670, 0.7814, 5.3338/), shape(c_2))
    test_num = 2
    call test_bulk_potential(c_2, a_2, expected_bulk_2, test_num)

    !Test 3 (Analytical) for bulk potential
    ! f(c) = c^2
    ! f'(c) = 2*c (expected_bulk)
    ! c = x^2 + y^2
    dx = 0.01
    dy = 0.01
    a_3 = (/0.0,0.0,1.0/)
    do i=1,3
        do j=1,3
            c_3(i,j) = ((i*dx)**2)+((j*dy)**2)
            expected_bulk_3(i,j) = 2*c_3(i,j)
        end do
    end do
    test_num = 3
    call test_bulk_potential(c_3,a_3,expected_bulk_3,test_num)
    !call test_bulk_potential(c_3,a_3,expeceted_bulk_3)

    !Test 1 for total potential
    dx = 1.0
    dy = 1.0
    kappa = 2.0
    expected_total_1 = reshape((/0.0, 12.0, 18.0, 8.0, 20.0, 62.0, -2.0, 10.0, 34.0/), shape(c_1))
    test_num = 1
    call test_total_potential(expected_bulk_1, c_1, dx, dx, kappa, expected_total_1, test_num,.false.)

    !Test 2 for total potential
    expected_total_2 = reshape((/-0.4902, 0.1380, 0.7517, 3.1837, 2.4189, 1.8470, 8.0138, -3.0764, -0.8056, 2.3637, &
                                 5.0620, 10.7708, -0.3965, 3.1670, -1.9585, 7.6538/), shape(c_2))
    test_num = 2
    call test_total_potential(expected_bulk_2, c_2, dx, dx, kappa, expected_total_2, test_num,.false.)

    !Test 3(analytical) for total potential
    dx = 0.01
    dy = 0.01
    kappa = 2.0
    test_num = 3
    expected_total_3 = expected_bulk_3 - 4.0*kappa
    call test_total_potential(expected_bulk_3,c_3,dx,dy,kappa,expected_total_3,test_num,.true.)
    !call test_total_potential(mu_3,c_3,dx,dx,kappa,expected_total_3)

    print*,'########################################################'
    print *, 'Starting free energy unit testing.'
    print*,'--------------------------------------------------------'


    !Test 1 for bulk free energy
    test_num = 1
    expected_fb_1 = reshape((/6.0, 17.0, 34.0, 17.0, 34.0, 121.0, 6.0, 17.0, 57.0/), shape(c_1))
    call test_bulk_free(c_1, a_1, expected_fb_1, test_num)

    !Test 2 for bulk free energy
    test_num = 2
    expected_fb_2 = reshape((/0.2945, 0.3524, 0.4534, 0.6789, 0.4958, 0.6164, 1.6139, 0.2671, 0.3619, 0.6789, 1.2783, &
                                2.3815, 0.3432, 0.6164, 0.3096, 1.6139/), shape(c_2))
    call test_bulk_free(c_2, a_2, expected_fb_2, test_num)

    !Test 1 for total free energy
    dx = 1.0
    dy = 1.0
    kappa = 0.34
    test_num = 1
    expected_total_F_1 = 311.71999

    call test_total_free_energy(c_1, expected_fb_1, dx, dy, kappa, expected_total_F_1, test_num)

    !Test 2 for total free energy
    dx = 0.01
    dy = 0.01
    kappa = 2.3
    test_num = 2
    expected_total_F_2 = 0.939175

    call test_total_free_energy(c_2, expected_fb_2, dx, dy, kappa, expected_total_F_2, test_num)

    print*,'########################################################'
    print*, "Starting Convergence test for F"
    print*,'--------------------------------------------------------'

    call free_energy_convergence(c_2, expected_fb_2, kappa)

    print*,'########################################################'
    print*, 'Testing First order derivatives - (for Chain rule term)'
    print*,'--------------------------------------------------------'

    Nx = 128
    Ny = 128

    call test_x_der(Nx, Ny)
    call test_y_der(Nx, Ny)



    print*,'########################################################'
    print*, "Starting testing for Cahn Hilliard solver"
    print*,'--------------------------------------------------------'

    Nx = 3
    Ny = 3

    dx = 0.01
    dy = 0.01
    dt = 1e-12 ! 1 picosecond timestep

    t_max = 1

    Kappa = 1.6

    M = 1.5

    C_mean = 0.7
    C_std = 0.1

    ! Try a 4th order polynomial
    allocate(a(6))
    a = (/5.0,2.1,2.2,2.3,2.4,2.5/)

    ! Test on del_Q subroutine

    ! Allocate dQ, Q_test, and Q_expected
    allocate(dQ(Nx,Ny))
    allocate(Q_test(Nx,Ny))
    allocate(dQ_expected(Nx,Ny))

    ! Set test values for Q_test
    do i = 1, Nx
        do j = 1, Ny
            Q_test(i,j) = ((i*dx)**2)+((j*dy)**2)
        end do
    end do

    call test_del_Q(dQ,Q_test,dQ_expected,dx,dy,Nx,Ny)

    ! Test on time_evolution subroutine

    ! Set dQ to zero
    dQ = 0.0

    ! Allocate c_test, c_expected, c_new, Q, and mu
    allocate(c_test(Nx,Ny))
    allocate(c_expected(Nx,Ny))
    allocate(c_new(Nx,Ny))

    ! initialize concentration grid with test values
    c_test = reshape((/0.82600226343913463,0.73763425858679454,0.75563253296536192, &
                       0.74015054492818821,0.73900628101493215,0.68999598441730403, &
                       0.70885112417583473,0.65750318317463130,0.80396934889945193/),shape(c_test))

    c_expected = reshape((/0.825592729545098,0.737684308382165,0.755776463302995, &
                           0.740158839048355,0.738720099275076,0.690382929167892, &
                           0.709132521563533,0.657924059121464,0.803700028541638/),shape(c_expected))

    call test_time_evolution(c_new,c_test,c_expected,Nx,Ny,dx,dy,dt,a,Kappa,M)
    

    print*,'########################################################'
    print*, "Starting grid test"
    print*,'--------------------------------------------------------'

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
    a = (/5.0,2.1,2.2,2.3,2.4,2.5/)

    ! Allocate grid
    allocate(c_grid(Nx,Ny))

    ! TESTS

    call test_grid_init(c_grid,Nx,Ny,C_mean,C_std,mean,std)

    call test_rand_normal(mean,std,C_mean,C_std)

    call test_stdnormal(mean,std)

    print*,'########################################################'
    print*, 'Tests Complete'
    

end program main
