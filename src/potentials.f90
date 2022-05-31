!> Potentials module contains routines to 
!! calculate bulk and total chemical potentials
module potentials

    use iso_fortran_env
    use grid

    implicit none

contains

    !> Subroutine to calculate bulk chemical potential.
    !! @param mu 2D grid to store bulk chemical potentials
    !! @param c 2D concentration grid
    !! @param a 1D array storing coefficients for bulk free energy
    subroutine bulk_potential(mu, c, T)

        real(real64) :: mu(:, :)
        real(real64), intent(in) :: c(:, :)
        real(real64), intent(in) :: T(:, :)
        real(real64) :: val
        integer :: nx, ny, n, i, j, k
        real(real64), parameter :: R = 8.314

        !Get mu grid size and number of user input coefficients
        nx = size(mu, 1)
        ny = size(mu, 2)

        !Initialize mu array with zero values
        mu = 0.0

        print*, c(1,6)

        ! Loop filling the mu array
        !$omp parallel do default(shared) private(j,i,k)
        do j = 1, ny
            do i = 1, nx
                mu(i, j) = mu(i, j) + R*T(i,j) * (log(c(i,j)) - log(1-c(i,j)) + 2)
                mu(i,j) =  mu(i,j) + 15911 + 3.335*T(i,j) - 2.0*(15911 + 3.335*T(i,j))*c(i,j)
            end do
        end do
        !$omp end parallel do

    end subroutine bulk_potential

    !> Subroutine to calculate total chemical potential.
    !! @param Q 2D grid to store total chemical potential
    !! @param mu 2D grid containing bulk chemical potential
    !! @param c 2D concentration grid
    !! @param dx Spatial grid spacing in x-direction
    !! @param dy Spatial grid spacing in y-direction
    !! @param Kappa Free energy gradient paramater
    !! @param conc_halo Concentration halo storing neighbour rank boundary data
    subroutine total_potential(Q, mu, c, T, dx, dy,conc_halo)

        real(real64), intent(in) :: mu(:,:)
        real(real64), intent(in) :: c(:,:)
        real(real64), intent(in) :: T(:,:)
        real(real64), intent(in) :: dx, dy
        real(real64), intent(in) :: conc_halo(:,:)
        real(real64) :: Q(:,:)
        real(real64) :: lap_x, lap_y,dx2,dy2
        integer :: nx, ny, i, j
        real(real64), parameter :: r_0 = 3.012e-10
        real(real64), dimension(:, :), allocatable ::  Kappa

        ! Get Q grid size
        nx = size(Q, 1)
        ny = size(Q, 2)

        allocate(Kappa(nx, ny))

         do j = 1, ny
            do i = 1, nx

                Kappa(i,j) = 0.5*r_0*r_0*c(i,j)*(1.0-c(i,j))*15911 + 3.335*T(i,j)
            
            end do
        end do


        ! Pre-compute 1.0/(dx^2) and 1.0/(dy^2)
        dx2 = 1.0 / (dx * dx)
        dy2 = 1.0 / (dy * dy)

        ! Start computing derivatives
        ! Corner nodes
        ! Top-left
        lap_x = (c(2, 1) - 2.0 * c(1, 1) + conc_halo(1,left)) * dx2
        lap_y = (c(1, 2) - 2.0 * c(1, 1) + conc_halo(1,up)) * dy2
        Q(1, 1) = mu(1, 1) - Kappa(1,1) * (lap_x + lap_y)

        ! Bottom-left
        lap_x = (c(2, ny) - 2.0 * c(1, ny) + conc_halo(ny, left)) * dx2
        lap_y = (conc_halo(1, down) - 2.0 * c(1, ny) + c(1, ny-1)) * dy2
        Q(1, ny) = mu(1, ny) - Kappa(1,ny) * (lap_x + lap_y)

        ! Bottom-right
        lap_x = (conc_halo(ny,right) - 2.0 * c(nx, ny) + c(nx-1, ny)) * dx2
        lap_y = (conc_halo(nx,down) - 2.0 * c(nx, ny) + c(nx, ny-1)) * dy2
        Q(nx, ny) = mu(nx, ny) - Kappa(nx,ny) * (lap_x + lap_y)

        ! Top-right
        !print*, conc_halo(nx,up)
        lap_x = (conc_halo(1, right) - 2.0 * c(nx,1 ) + c(nx-1,1 )) * dx2
        lap_y = (c(nx, 2) - 2.0 * c(nx, 1) + conc_halo(nx,up)) * dy2
        Q(nx, 1) = mu(nx, 1) - Kappa(nx,1) * (lap_x + lap_y)

        ! Boundary nodes
        ! Top - j = 1
        !$omp parallel do default(shared) private(i,lap_x,lap_y)
        do i = 2, nx - 1
            lap_y = (c(i, 2) - 2.0 * c(i, 1) + conc_halo(i, up)) * dy2
            lap_x = (c(i + 1, 1) - 2.0 * c(i, 1) + c(i - 1, 1)) * dx2
            Q(i, 1) = mu(i, 1) - Kappa(i,1) * (lap_x + lap_y)
            lap_y = (conc_halo(i, down) - 2.0 * c(i, ny) + c(i, ny - 1)) * dy2
            lap_x = (c(i + 1, ny) - 2.0 * c(i, ny) + c(i - 1, ny)) * dx2
            Q(i, ny) = mu(i, ny) - Kappa(i,ny) * (lap_x + lap_y)
        end do
        !$omp end parallel do
        ! LHS - i = 1
        !$omp parallel do default(shared) private(j,lap_x,lap_y)
        do j = 2, ny - 1
            lap_y = (c(1, j + 1) - 2.0 * c(1, j) + c(1, j - 1)) * dy2
            lap_x = (c(2, j) - 2.0 * c(1, j) + conc_halo(j, left)) * dx2
            Q(1, j) = mu(1,j) - Kappa(1,j) * (lap_x + lap_y)
            lap_y = (c(nx, j + 1) - 2.0 * c(nx, j) + c(nx, j - 1)) * dy2
            lap_x = (conc_halo(j,right) - 2.0 * c(nx, j) + c(nx - 1, j)) * dx2
            Q(nx, j) = mu(nx, j) - Kappa(nx,j) * (lap_x + lap_y)
        end do
        !$omp end parallel do
        ! Bulk (non-boundary) nodes
        !$omp parallel do default(shared) private(j,i,lap_x,lap_y)
        do i = 2, nx - 1
            do j = 2, ny - 1
                lap_y = (c(i, j + 1) - 2 * c(i, j) + c(i, j - 1)) * dy2
                lap_x = (c(i + 1, j) - 2 * c(i, j) + c(i - 1, j)) * dx2
                Q(i, j) = mu(i, j) - Kappa(i,j) * (lap_x + lap_y)
            end do
        end do
        !$omp end parallel do

    end subroutine total_potential

end module potentials
