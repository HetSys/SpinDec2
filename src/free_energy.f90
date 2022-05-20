module free_energy

    use iso_fortran_env
    use grid

    implicit none

contains

    !******************************************************************************
    !> bulk_free_energy
    !!
    !! subroutine to calculate Bulk Free Energy at time t
    !!
    !!@param f_b: bulk free energy array at time t
    !!@param c: 2D concentration grid - c(x,y)
    !!@param a: 1D array storing coefficients provided by user
    !!@param kappa : gradient term coefficient
    !*****************************************************************************

    subroutine bulk_free_energy(f_b, c, a)

        real(real64), intent(out), allocatable :: f_b(:, :)
        real(real64), intent(in) :: c(:, :)
        real(real64), intent(in) :: a(:)
        integer :: nx, ny, n
        integer :: i, j, k

        !Get f_b grid size and number of user input coefficients
        Nx = size(c, 1)
        Ny = size(c, 2)
        n = size(a)

        f_b = 0

        !Loop filling the f_b array
        !$omp parallel do default(shared) private(j,i,k)
        do j = 1, Ny
            do i = 1, Nx
                do k = 1, n
                    f_b(i, j) = f_b(i, j) + a(k) * c(i, j)**(k - 1)
                end do
            end do
        end do

    end subroutine

    !******************************************************************************
    !> total_free_energy
    !!
    !! subroutine to calculate Total Free Energy at time t
    !!
    !!@param F: total free energy at time t
    !!@param c: 2D concentration grid - c(x,y)
    !!@param dx,dy: spatial step sizes in x and y
    !!@param kappa : gradient term coefficient
    !*****************************************************************************

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
        grad_x = (c(1, 2) -  conc_halo(1, left)) * dx2
        grad_y = (c(2, 1) -  conc_halo(1, up)) * dy2
        P = grad_x*grad_x + grad_y*grad_y
        F = F + (f_b(1, 1) + 0.5 * kappa * P ) * dx * dy

        ! Bottom-left
        grad_x = (c(nx, 2) - conc_halo(nx,left)) * dx2
        grad_y = (conc_halo(1, down) - c(nx - 1, 1)) * dy2
        P = grad_x*grad_x + grad_y*grad_y
        F = F + (f_b(nx, 1) + 0.5 * kappa * P ) * dx * dy

        ! Bottom-right
        grad_x = (conc_halo(nx, right) - c(nx, ny - 1)) * dx2
        grad_y = (conc_halo(nx, down) - c(nx - 1, ny)) * dy2
        P = grad_x*grad_x + grad_y*grad_y
        F = F + (f_b(nx, ny) + 0.5 * kappa * P ) * dx * dy

        ! Top-right
        grad_x = (conc_halo(1, right) - c(1, ny - 1)) * dx2
        grad_y = (c(2, ny) - conc_halo(nx, up)) * dy2
        P = grad_x*grad_x + grad_y*grad_y
        F = F + (f_b(1, ny) + 0.5 * kappa * P ) * dx * dy

        ! Boundary nodes
        ! LHS - j = 1
        do i = 2, nx - 1
            grad_x = (c(i, 2) - conc_halo(i,left)) * dx2
            grad_y = (c(i + 1, 1) - c(i - 1, 1)) * dy2
            P = grad_x*grad_x + grad_y*grad_y
            F = F + (f_b(i, 1) + 0.5 * kappa * P ) * dx * dy
        end do

        ! RHS - j = ny
        do i = 2, nx - 1
            grad_x = (conc_halo(i, right) - c(i, ny - 1)) * dx2
            grad_y = (c(i + 1, ny) - c(i - 1, ny)) * dy2
            P = grad_x*grad_x + grad_y*grad_y
            F = F + (f_b(i, ny) + 0.5 * kappa * P ) * dx * dy
        end do

        ! Top - i = 1
        do j = 2, ny - 1
            grad_x = (c(1, j + 1) - c(1, j - 1)) * dx2
            grad_y = (c(2, j) - conc_halo(j,up)) * dy2
            P = grad_x*grad_x + grad_y*grad_y
            F = F + (f_b(1, j) + 0.5 * kappa * P ) * dx * dy
        end do

        ! Bottom - i = nx
        do j = 2, ny - 1
            grad_x = (c(nx, j + 1) - c(nx, j - 1)) * dx2
            grad_y = (conc_halo(j,down) - c(nx - 1, j)) * dy2
            P = grad_x*grad_x + grad_y*grad_y
            F = F + (f_b(nx, j) + 0.5 * kappa * P ) * dx * dy
        end do

        ! Bulk (non-boundary) nodes
        !$omp parallel do default(shared) private(j,i,grad_x,grad_y,P) reduction(+:F)
        do j = 2, ny - 1
            do i = 2, nx - 1
                grad_x = (c(i, j + 1) -  c(i, j - 1)) * dx2
                grad_y = (c(i + 1, j) -  c(i - 1, j)) * dy2
                P = grad_x*grad_x + grad_y*grad_y
                F = F + (f_b(i, j) + 0.5 * kappa * P ) * dx * dy
            end do
        end do

    end subroutine

end module free_energy
