module potentials

    use iso_fortran_env

    implicit none

contains

    ! Subroutine to calculate bulk chemical potential
    ! @param mu: 2D grid to store bulk chemical potentials
    ! @param c: 2D concentration grid
    ! @param a: 1D array storing coefficients provided by user
    subroutine bulk_potential(mu, c, a)

        real(real64) :: mu(:, :)
        real(real64), intent(in) :: c(:, :)
        real(real64), intent(in) :: a(:)
        integer :: nx, ny, n, i, j, k

        !Get mu grid size and number of user input coefficients
        nx = size(mu, 1)
        ny = size(mu, 2)
        n = size(a)

        !Initialize mu array with zero values
        mu = 0.0

        ! Loop filling the mu array
        !$omp parallel do default(shared) private(j,i,k)
        do j = 1, ny
            do i = 1, nx
                do k = 1, n - 1
                    mu(i, j) = mu(i, j) + k * a(k + 1) * c(i, j)**(k - 1)
                end do
            end do
        end do

    end subroutine bulk_potential

    ! Subroutine to calculate total chemical potential
    ! @param Q: 2D grid to store total chemical potential
    ! @param mu: 2D grid containing bulk chemical potential
    ! @param c: 2D concentration grid
    ! @param dx: spatial grid spacing in x-direction
    ! @param dy: spatial grid spacing in y-direction
    ! @param Kappa: Free energy gradient paramater
    subroutine total_potential(Q, mu, c, dx, dy, Kappa,conc_halo)

        real(real64), intent(in) :: mu(:,:)
        real(real64), intent(in) :: c(:,:)
        real(real64), intent(in) :: dx, dy, Kappa
        real(real64), intent(in) :: conc_halo(:,:)
        real(real64) :: Q(:,:)
        real(real64) :: lap_x, lap_y,dx2,dy2
        integer :: nx, ny, i, j

        ! Get Q grid size
        nx = size(Q, 1)
        ny = size(Q, 2)
        

        ! Pre-compute 1.0/(dx^2) and 1.0/(dy^2)
        dx2 = 1.0 / (dx * dx)
        dy2 = 1.0 / (dy * dy)

        ! Start computing derivatives
        ! Corner nodes
        ! Top-left
        lap_x = (c(2, 1) - 2.0 * c(1, 1) + c(nx, 1)) * dx2
        lap_y = (c(1, 2) - 2.0 * c(1, 1) + c(1, ny)) * dy2
        Q(1, 1) = mu(1, 1) - Kappa * (lap_x + lap_y)

        ! Bottom-left
        lap_x = (c(1, 1) - 2.0 * c(nx, 1) + c(nx - 1, 1)) * dx2
        lap_y = (c(nx, 2) - 2.0 * c(nx, 1) + c(nx, ny)) * dy2
        Q(nx, 1) = mu(nx, 1) - Kappa * (lap_x + lap_y)

        ! Bottom-right
        lap_x = (c(1, ny) - 2.0 * c(nx, ny) + c(nx - 1, ny)) * dx2
        lap_y = (c(nx, 1) - 2.0 * c(nx, ny) + c(nx, ny - 1)) * dy2
        Q(nx, ny) = mu(nx, ny) - Kappa * (lap_x + lap_y)

        ! Top-right
        lap_x = (c(2, ny) - 2.0 * c(1, ny) + c(nx, ny)) * dx2
        lap_y = (c(1, 1) - 2.0 * c(1, ny) + c(1, ny - 1)) * dy2
        Q(1, ny) = mu(1, ny) - Kappa * (lap_x + lap_y)

        ! Boundary nodes
        ! LHS - j = 1
        do i = 2, nx - 1
            lap_x = (c(i + 1, 1) - 2.0 * c(i, 1) + c(i - 1, 1)) * dx2
            lap_y = (c(i, 2) - 2.0 * c(i, 1) + c(i, ny)) * dy2
            Q(i, 1) = mu(i, 1) - Kappa * (lap_x + lap_y)
        end do

        ! RHS - j = ny
        do i = 2, nx - 1
            lap_x = (c(i + 1, ny) - 2.0 * c(i, ny) + c(i - 1, ny)) * dx2
            lap_y = (c(i, 1) - 2.0 * c(i, ny) + c(i, ny - 1)) * dy2
            Q(i, ny) = mu(i, ny) - Kappa * (lap_x + lap_y)
        end do

        ! Top - i = 1
        do j = 2, ny - 1
            lap_x = (c(2, j) - 2.0 * c(1, j) + c(nx, j)) * dx2
            lap_y = (c(1, j + 1) - 2.0 * c(1, j) + c(1, j - 1)) * dy2
            Q(1, j) = mu(1,j) - Kappa * (lap_x + lap_y)
        end do

        ! Bottom - i = nx
        do j = 2, ny - 1
            lap_x = (c(1, j) - 2.0 * c(nx, j) + c(nx - 1, j)) * dx2
            lap_y = (c(nx, j + 1) - 2.0 * c(nx, j) + c(nx, j - 1)) * dy2
            Q(nx, j) = mu(nx, j) - Kappa * (lap_x + lap_y)
        end do

        ! Bulk (non-boundary) nodes
        !$omp parallel do default(shared) private(j,i,lap_x,lap_y)
        do i = 2, nx - 1
            do j = 2, ny - 1
                lap_x = (c(i + 1, j) - 2 * c(i, j) + c(i - 1, j)) * dx2
                lap_y = (c(i, j + 1) - 2 * c(i, j) + c(i, j - 1)) * dy2
                Q(i, j) = mu(i, j) - Kappa * (lap_x + lap_y)
            end do
        end do

        ! !Loop filling the Q array
        ! do j = 0, ny - 1
        !     do i = 0, nx - 1
        !         !X and Y double derivatives using finite differences
        !         !Peridic boundary conditions are taken care using mod function
        !         x_lap = (c(modulo(i + 1, nx), j) - 2 * c(i, j) + c(modulo(i - 1, nx), j)) / (dx * dx)
        !         y_lap = (c(i, modulo(j + 1, ny)) - 2 * c(i, j) + c(i, modulo(j - 1, ny))) / (dy * dy)
        !         Q(i, j) = mu(i, j) - Kappa * (x_lap + y_lap)
        !     end do
        ! end do

    end subroutine total_potential

end module potentials
