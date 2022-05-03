program main
    
    use iso_fortran_env
    use cahn_hilliard
    use grid_modules
    use potentials
    
    implicit none

    real(real64), dimension(:,:), allocatable :: grid ! conc. grid
    real(real64), dimension(:,:), allocatable :: c_new ! new conc. grid
    real(real64), dimension(:,:), allocatable :: mu   ! bulk chem. pot.
    real(real64), dimension(:,:), allocatable :: Q    ! total chem. pot.
    real(real64), dimension(:,:), allocatable :: dQ   ! 2nd derivative of Q
    real(real64), dimension(:), allocatable :: a
    real(real64) :: C_mean, C_std
    real(real64) :: dx, dy, dt ! spatial and temporal grid spacings
    real(real64) :: Kappa ! free energy gradient parameter
    real(real64) :: M ! Mobility
    integer :: Nx, Ny
    integer :: t_max ! number of timesteps
    integer :: i,j,t ! counters
    character(len=100) :: out_file


    Nx = 50
    Ny = 50

    dx = 0.01
    dy = 0.01
    dt = 1e-12 ! 1 picosecond timestep

    t_max = 10000

    Kappa = 1.6

    M = 1.5

    C_mean = 0.7
    C_std = 0.1

    ! Try a 4th order polynomial
    allocate(a(6))

    a = (/5.0,2.1,2.2,2.3,2.4,2.5/)

    ! Allocate grid
    allocate(grid(Nx,Ny))

    ! Initialize grid
    call grid_init(grid,Nx,Ny,C_mean,C_std)

    ! Write original grid to file
    open (99, file = "grid.dat", status = "new")
    do i = 1, Nx
        write(99,*) grid(i,:)
    end do
    close(99)

    ! Allocate mu
    allocate(mu(Nx,Ny))
    ! Allocate Q
    allocate(Q(Nx,Ny))
    ! Allocate dQ
    allocate(dQ(Nx,Ny))
    ! Allocate c_new
    allocate(c_new(Nx,Ny))

    
    ! Grid evolution
    do t = 1, t_max
        ! Get total chemical potentials
        call total_potential(Q,mu,grid,dx,dy,Kappa)

        ! Get bulk chemical potentials
        call bulk_potential(mu,grid,a)

        ! Get second derivative of Q
        call del_Q(dQ,Q,M,dx,dy)

        ! Get new concentrations for current timesteps
        call time_evolution(grid,c_new,dQ,dt,Nx,Ny)

        ! set grid to c_new
        grid = c_new

        ! ! write new concentration grid to data file
        ! write(out_file, "(a,i0,a)") "grid_t",t,".dat"
        ! open (99, file = out_file, status = "new")
        ! do i = 1, Nx
        !     write(99,*) c_new(i,:)
        ! end do
        ! close(99)
    end do

        ! write final grid to data file
        !write(out_file, "(a,i0,a)") "grid_t",t,".dat"
        open (99, file = "grid_fin.dat", status = "new")
        do i = 1, Nx
            write(99,*) c_new(i,:)
        end do
        close(99)

end program main