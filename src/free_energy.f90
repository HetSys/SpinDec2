!> Free energy module contains routines to 
!! calculate bulk and total free energy
module free_energy

    use iso_fortran_env
    use grid

    implicit none

contains

    !> Subroutine to calculate Bulk Free Energy.
    !! @param f_b 2D grid for bulk free energy 
    !! @param c 2D concentration grid - c(x,y)
    !! @param a 1D array storing coefficients provided by user
    subroutine bulk_free_energy(f_b, c, a)

        real(real64), intent(out) :: f_b(:, :)
        real(real64), intent(in) :: c(:, :)
        real(real64), intent(in) :: a(:)
        integer :: nx, ny, n
        integer :: i, j, k
        real(real64) :: val

        !Get f_b grid size and number of user input coefficients
        Nx = size(c, 1)
        Ny = size(c, 2)
        n = size(a)

        f_b = 0

        !Loop filling the f_b array
        !$omp parallel do default(shared) private(j,i,k)
        do j = 1, Ny
            do i = 1, Nx
                val = 1.0_real64
                do k = 1, n
                    f_b(i, j) = f_b(i, j) + a(k) * val
                    val = val * c(i,j)
                end do
            end do
        end do
        !$omp end parallel do

    end subroutine bulk_free_energy

    !> Subroutine to calculate Total Free Energy
    !! @param F total free energy at time t
    !! @param c 2D concentration grid - c(x,y)
    !! @param f_b 2D grid storing bulk free energy 
    !! @param dx Spatial step size in x
    !! @param dy Spatial step size in y
    !! @param kappa Free energy gradient paramater
    !! @param conc_halo Concentration halo storing neighbour rank boundary data
    subroutine total_free_energy(F, c, f_b, dx, dy, kappa,conc_halo)

        real(real64), intent(in) :: c(:,:)
        real(real64), intent(in) :: conc_halo(:,:)
        real(real64) :: f_b(:,:)
        integer :: nx, ny
        integer :: i, j
        real(real64), intent(out) :: F
        real(real64) :: dx, dy, dx2, dy2
        real(real64) :: kappa
        real(real64) :: P, grad_x, grad_y

        !Get mu grid size and number of user input coefficients
        nx = size(c, 1)
        ny = size(c, 2)

        F = 0.0
        ! Pre-compute 1.0/(dx^2) and 1.0/(dy^2)
        dx2 = 1.0 / (2.0 * dx)
        dy2 = 1.0 / (2.0 * dy)

        ! Start computing derivatives
        ! Corner nodes
        ! Top-left
        grad_x = (c(2, 1) -  conc_halo(1, left)) * dx2
        grad_y = (c(1, 2) -  conc_halo(1, up)) * dy2
        P = grad_x*grad_x + grad_y*grad_y
        F = F + (f_b(1, 1) + 0.5 * kappa * P ) * dx * dy

        ! Bottom-left
        grad_x = (c(2, ny) - conc_halo(ny,left)) * dx2
        grad_y = (conc_halo(1, down) - c(1,ny - 1)) * dy2
        P = grad_x*grad_x + grad_y*grad_y
        F = F + (f_b(1, ny) + 0.5 * kappa * P ) * dx * dy

        ! Bottom-right
        grad_x = (conc_halo(ny, right) - c(nx-1, ny)) * dx2
        grad_y = (conc_halo(nx, down) - c(nx, ny-1)) * dy2
        P = grad_x*grad_x + grad_y*grad_y
        F = F + (f_b(nx, ny) + 0.5 * kappa * P ) * dx * dy

        ! Top-right
        grad_x = (conc_halo(1, right) - c(nx - 1,1)) * dx2
        grad_y = (c(nx,2) - conc_halo(nx, up)) * dy2
        P = grad_x*grad_x + grad_y*grad_y
        F = F + (f_b(nx, 1) + 0.5 * kappa * P ) * dx * dy

        ! Boundary nodes
        ! Top - j = 1
        !$omp parallel do default(shared) private(i,grad_x,grad_y,P) reduction(+:F)
        do i = 2, nx - 1
            grad_y = (c(i, 2) - conc_halo(i,up)) * dy2
            grad_x = (c(i + 1, 1) - c(i - 1, 1)) * dx2
            P = grad_x*grad_x + grad_y*grad_y
            F = F + (f_b(i, 1) + 0.5 * kappa * P ) * dx * dy
            grad_y = (conc_halo(i, down) - c(i, ny - 1)) * dy2
            grad_x = (c(i + 1, ny) - c(i - 1, ny)) * dx2
            P = grad_x*grad_x + grad_y*grad_y
            F = F + (f_b(i, ny) + 0.5 * kappa * P ) * dx * dy
        end do
        !$omp end parallel do
        ! Left - i = 1
        !$omp parallel do default(shared) private(j,grad_x,grad_y,P) reduction(+:F)
        do j = 2, ny - 1
            grad_y = (c(1, j + 1) - c(1, j - 1)) * dy2
            grad_x = (c(2, j) - conc_halo(j,left)) * dx2
            P = grad_x*grad_x + grad_y*grad_y
            F = F + (f_b(1, j) + 0.5 * kappa * P ) * dx * dy
            grad_y = (c(nx, j + 1) - c(nx, j - 1)) * dy2
            grad_x = (conc_halo(j,right) - c(nx - 1, j)) * dx2
            P = grad_x*grad_x + grad_y*grad_y
            F = F + (f_b(nx, j) + 0.5 * kappa * P ) * dx * dy
        end do
        !$omp end parallel do
        ! Bulk (non-boundary) nodes
        !$omp parallel do default(shared) private(j,i,grad_x,grad_y,P) reduction(+:F)
        do j = 2, ny - 1
            do i = 2, nx - 1
                grad_y = (c(i, j + 1) -  c(i, j - 1)) * dy2
                grad_x = (c(i + 1, j) -  c(i - 1, j)) * dx2
                P = grad_x*grad_x + grad_y*grad_y
                F = F + (f_b(i, j) + 0.5 * kappa * P ) * dx * dy
            end do
        end do
        !$omp end parallel do
    end subroutine total_free_energy

end module free_energy