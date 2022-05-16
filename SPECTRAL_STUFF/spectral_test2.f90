!gfortran -std=f2018 -fbounds-check src/Input_parser.f90 src/checkpoint.f90 src/free_energy.f90 src/grid_mod.f90 src/potentials.f90 src/io.f90 src/cahn_hilliard.f90 src/main.f90 `nf-config --fflags --flibs`

program main
    use spectral
    use iso_fortran_env
    use grid
    use potentials
    use io


    implicit none

    real(kind=real64) ,dimension(:,:), allocatable:: c_in
    real(kind=real64) ,dimension(:), allocatable:: f_tot
    real(kind=real64) ,dimension(:,:,:), allocatable:: c
    real(kind=real64) ,dimension(:,:), allocatable:: c_in_p
    real(kind=real64) ,dimension(:,:), allocatable:: c_out
    real(kind=real64) ,dimension(:), allocatable:: coeffs
    real(kind=real64) :: M,kappa,dt,c0,c_std,ma,mb
    integer :: k,nx,ny,nt,nc





    nc = 5


    ma = 1
    mb = 1


    kappa = 1
    dt = 0.02
    c0 = 0.7
    c_std = 0.3

    M = (MB * (1 - c0) + MA/MB * c0) * c0 * (1 - c0)
    nx = 50
    ny = 50
    nt = 100000

    allocate(c(nx,ny,nt))
    allocate(c_in(3,3))
    allocate(c_in_p(3,3))
    allocate(c_out(nx,ny))
    allocate(coeffs(5))
    allocate(f_tot(nt))

    coeffs(1) = 100
    coeffs(2) = 0
    coeffs(3) = 1.0*1.7
    coeffs(4) = -2.0*1.7
    coeffs(5) = 1.0*1.7

    f_tot = 0

    call get_seed(134254)
    ! Initialize grid
    call grid_init(c(:, :, 1), Nx, Ny, c0, c_std)
    !c(:,:,1)=0
    !c(10:,10:,1)=1
    call spectral_method_iter(c(:, :, 1),c(:, :, 1),coeffs,dt,M,Kappa,c_out,1)

    c(:,:,2) = c_out
    ! Grid evolution
    do k = 2, Nt-1
        call spectral_method_iter(c(:, :, k),c(:, :, k-1),coeffs,dt,M,Kappa,c_out,0)
        c(:,:,k+1)=c_out
    end do

    !Writer for using constant M
    call write_netcdf(c, F_tot, coeffs, Nc, Nx, Ny, Nt, dt, c0, MA, MB, kappa)



    !print*,c_out
end program main
