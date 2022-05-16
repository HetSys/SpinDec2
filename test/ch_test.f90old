module ch_test

    use iso_fortran_env
    use grid
    use cahn_hilliard
    use potentials

    implicit none

    contains

    subroutine test_del_Q(dQ,Q_test,dQ_expected,M,dx,dy,Nx,Ny)

        integer, intent(in) :: Nx, Ny
        real(real64), dimension(Nx,Ny), intent(in) :: Q_test, dQ_expected
        real(real64), dimension(Nx,Ny), intent(out) :: dQ
        real(real64), intent(in) :: M, dx, dy
        logical :: result = .TRUE.
        integer :: i,j ! counters

        call del_Q(dQ,Q_test,M,dx,dy,Nx,Ny)

        do i = 1, Nx
            do j = 1 , Ny
                ! print*, abs(dQ(i,j)-dQ_expected(i,j))
                if(abs(dQ(i,j)-dQ_expected(i,j)) > 1e-2) then
                    result = .FALSE.
                    exit ! do not continue checking
                end if
            end do
        end do
        
        if (result) then
            print*, "del_Q passes test"
        else
            print*, "del_Q does not pass test"
        end if

    end subroutine test_del_Q


    subroutine test_time_evolution(c_new,c_test,c_expected,Nx,Ny,dx,dy,dt,a,Kappa,M)

        integer, intent(in) :: Nx, Ny
        real(real64), dimension(Nx,Ny), intent(in) :: c_test, c_expected
        real(real64), intent(in) :: dx, dy, dt, Kappa, M
        real(real64), dimension(:), intent(in) :: a
        real(real64), dimension(Nx,Ny), intent(out) :: c_new
        real(real64), dimension(Nx,Ny) :: mu, Q, dQ
        logical :: result = .TRUE.
        integer :: i,j ! counters


        ! Get bulk chemical potential
        call bulk_potential(mu,c_test,a)
        ! Get total chemical potential
        call total_potential(Q,mu,c_test,dx,dy,Kappa)
        ! Get dQ
        call del_Q(dQ,Q,M,dx,dy,Nx,Ny)
        ! Get grid at next time step
        call time_evolution(c_test,c_new,dQ,dt,Nx,Ny)

        call del_Q(dQ,Q,M,dx,dy,Nx,Ny)

        do i = 1, Nx
            do j = 1 , Ny
                ! print *, abs(c_new(i,j)-c_expected(i,j))
                if(abs(c_new(i,j)-c_expected(i,j)) > 1e-3) then
                    result = .FALSE.
                    exit ! do not continue checking
                end if
            end do
        end do
        
        if (result) then
            print*, "time_evolution passes test"
        else
            print*, "time_evolution does not pass test"
        end if

    end subroutine test_time_evolution

end module ch_test

program main

    use iso_fortran_env
    use ch_test

    implicit none

    real(real64), dimension(:,:), allocatable :: c_new ! new conc. grid
    real(real64), dimension(:,:), allocatable :: dQ   ! 2nd derivative of Q
    real(real64), dimension(:), allocatable :: a
    real(real64) :: C_mean, C_std
    real(real64) :: dx, dy, dt ! spatial and temporal grid spacings
    real(real64) :: Kappa ! free energy gradient parameter
    real(real64) :: M ! Mobility
    integer :: Nx, Ny
    integer :: t_max ! number of timesteps

    real(real64), dimension(:,:), allocatable :: Q_test, dQ_expected
    real(real64), dimension(:,:), allocatable :: c_test, c_expected

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

    ! Set test values for Q
    Q_test = reshape((/0.905,0.553,0.633, &
                       0.788,0.770,0.781, &
                       0.735,0.672,0.653/),shape(Q_test))
    ! Q_test = reshape((/0.826002263439135,0.740150544928188,0.708851124175835, &
    !                    0.737634258586795,0.739006281014932,0.657503183174631, &
    !                    0.755632532965362,0.689995984417304,0.803969348899452/),shape(Q_test))

    ! Set known del_Q values
    dQ_expected = reshape((/-13665.0,11520.0,5400.0, &
                             585.0,-4290.0,-4200.0, &
                             1170.0,345.0,3135.0/),shape(dQ_expected))

    ! dQ_expected = reshape((/-5426.10889650539,48.8021000167921,2883.31265607101, &
    !                          414.108393703222,-1961.11729429216,5391.72419967732, &
    !                          526.075852218563,4181.62155208077,-6058.41856297014/),shape(dQ_expected))

    call test_del_Q(dQ,Q_test,dQ_expected,M,dx,dy,Nx,Ny)

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

end program main
