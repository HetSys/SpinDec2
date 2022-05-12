module cahn_hilliard

    use iso_fortran_env

    implicit none

contains

    subroutine del_Q(dQ, Q, dx, dy, Nx, Ny)
        ! Subroutine to perform spatial dervative of Q
        ! dQ is a 2D grid to store the derivatives of Q
        ! Q is the 2D grid of total chemical potentials
        ! M is the mobility of the species
        ! dx and dy are the spatial grid spacings in x and y


        integer, intent(in) :: Nx, Ny
        real(real64), dimension(Nx, Ny), intent(in) :: Q
        real(real64), dimension(Nx, Ny), intent(out) :: dQ
        real(real64), intent(in) :: dx, dy
        real(real64) :: dx2, dy2
        real(real64) :: der_x, der_y
        integer :: i, j ! counters

        ! Pre-compute 1.0/(dx^2) and 1.0/(dy^2)
        dx2 = 1.0 / (dx * dx)
        dy2 = 1.0 / (dy * dy)

        ! Start computing derivatives
        ! Corner nodes
        ! Top-left
        der_x = (Q(2, 1) - 2.0 * Q(1, 1) + Q(Nx, 1)) * dx2
        der_y = (Q(1, 2) - 2.0 * Q(1, 1) + Q(1, Ny)) * dy2
        dQ(1, 1) = (der_x + der_y)

        ! Bottom-left
        der_x = (Q(1, 1) - 2.0 * Q(Nx, 1) + Q(Nx - 1, 1)) * dx2
        der_y = (Q(Nx, 2) - 2.0 * Q(Nx, 1) + Q(Nx, Ny)) * dy2
        dQ(Nx, 1) = (der_x + der_y)

        ! Bottom-right
        der_x = (Q(1, Ny) - 2.0 * Q(Nx, Ny) + Q(Nx - 1, Ny)) * dx2
        der_y = (Q(Nx, 1) - 2.0 * Q(Nx, Ny) + Q(Nx, Ny - 1)) * dy2
        dQ(Nx, Ny) = (der_x + der_y)

        ! Top-right
        der_x = (Q(2, Ny) - 2.0 * Q(1, Ny) + Q(Nx, Ny)) * dx2
        der_y = (Q(1, 1) - 2.0 * Q(1, Ny) + Q(1, Ny - 1)) * dy2
        dQ(1, Ny) = (der_x + der_y)

        ! Boundary nodes
        ! LHS - j = 1
        do i = 2, Nx - 1
            der_x = (Q(i + 1, 1) - 2.0 * Q(i, 1) + Q(i - 1, 1)) * dx2
            der_y = (Q(i, 2) - 2.0 * Q(i, 1) + Q(i, Ny)) * dy2
            dQ(i, 1) = (der_x + der_y)
        end do

        ! RHS - j = Ny
        do i = 2, Nx - 1
            der_x = (Q(i + 1, Ny) - 2.0 * Q(i, Ny) + Q(i - 1, Ny)) * dx2
            der_y = (Q(i, 1) - 2.0 * Q(i, Ny) + Q(i, Ny - 1)) * dy2
            dQ(i, Ny) = (der_x + der_y)
        end do

        ! Top - i = 1
        do j = 2, Ny - 1
            der_x = (Q(2, j) - 2.0 * Q(1, j) + Q(Nx, j)) * dx2
            der_y = (Q(1, j + 1) - 2.0 * Q(1, j) + Q(1, j - 1)) * dy2
            dQ(1, j) = (der_x + der_y)
        end do

        ! Bottom - i = Nx
        do j = 2, Ny - 1
            der_x = (Q(1, j) - 2.0 * Q(Nx, j) + Q(Nx - 1, j)) * dx2
            der_y = (Q(Nx, j + 1) - 2.0 * Q(Nx, j) + Q(Nx, j - 1)) * dy2
            dQ(Nx, j) = (der_x + der_y)
        end do

        ! Bulk (non-boundary) nodes
        do i = 2, Nx - 1
            do j = 2, Ny - 1
                der_x = (Q(i + 1, j) - 2 * Q(i, j) + Q(i - 1, j)) * dx2
                der_y = (Q(i, j + 1) - 2 * Q(i, j) + Q(i, j - 1)) * dy2
                dQ(i, j) = (der_x + der_y)
            end do
        end do

    end subroutine del_Q


    subroutine dQ_dy(dQy,Q,dy,Nx,Ny)
        ! Subroutine to calculate dQ/dx
        ! Using central differences
        integer, intent(in) :: Nx, Ny
        real(real64), dimension(Nx,Ny), intent(in) :: Q
        real(real64), intent(in) :: dy
        real(real64), dimension(Nx,Ny), intent(out) :: dQy
        real(real64) :: dy_inv
        integer :: i,j ! counters

        dy_inv = 1.0/(2.0*dy)

        ! LHS and RHS Boundary
        do i = 1, Nx
            dQy(i,1) = (Q(i,2) - Q(i,Ny))*dy_inv    ! LHS
            dQy(i,Ny) = (Q(i,1) - Q(i,Ny-1))*dy_inv ! RHS
        end do

        ! Bulk (non-boundary nodes)
        do j = 2, Ny-1
            do i = 1, Nx
                dQy(i,j) = (Q(i,j+1) - Q(i,j-1))*dy_inv
            end do
        end do

    end subroutine dQ_dy

    subroutine dQ_dx(dQx,Q,dx,Nx,Ny)
        ! Subroutine to calculate dQ/dx
        ! Using central differences
        integer, intent(in) :: Nx, Ny
        real(real64), dimension(Nx,Ny), intent(in) :: Q
        real(real64), intent(in) :: dx
        real(real64), dimension(Nx,Ny), intent(out) :: dQx
        real(real64) :: dx_inv
        integer :: i,j ! counters

        dx_inv = 1.0/(2.0*dx)

        ! Top and Bottom Boundary
        do j = 1, Ny
            dQx(1,j) = (Q(2,j) - Q(Nx,j))*dx_inv    ! Top
            dQx(Nx,j) = (Q(1,j) - Q(Nx-1,j))*dx_inv ! Bottom
        end do 

        ! Bulk (non-boundary nodes)
        do j = 1, Ny
            do i = 2, Nx-1
                dQx(i,j) = (Q(i+1,j) - Q(i-1,j))*dx_inv
            end do
        end do

    end subroutine dQ_dx


    subroutine dM_dy(dMy,M,dy,Nx,Ny)
        ! Subroutine to calculate dM/dx
        ! Using central differences
        integer, intent(in) :: Nx, Ny
        real(real64), dimension(Nx,Ny), intent(in) :: M
        real(real64), intent(in) :: dy
        real(real64), dimension(Nx,Ny), intent(out) :: dMy
        real(real64) :: dy_inv
        integer :: i,j ! counters

        dy_inv = 1.0/(2.0*dy)

        ! LHS and RHS Boundary
        do i = 1, Nx
            dMy(i,1) = (M(i,2) - M(i,Ny))*dy_inv    ! LHS
            dMy(i,Ny) = (M(i,1) - M(i,Ny-1))*dy_inv ! RHS
        end do

        ! Bulk (non-boundary nodes)
        do j = 2, Ny-1
            do i = 1, Nx
                dMy(i,j) = (M(i,j+1) - M(i,j-1))*dy_inv
            end do
        end do

    end subroutine dM_dy

    subroutine dM_dx(dMx,M,dx,Nx,Ny)
        ! Subroutine to calculate dM/dx
        ! Using central differences
        integer, intent(in) :: Nx, Ny
        real(real64), dimension(Nx,Ny), intent(in) :: M
        real(real64), intent(in) :: dx
        real(real64), dimension(Nx,Ny), intent(out) :: dMx
        real(real64) :: dx_inv
        integer :: i,j ! counters

        dx_inv = 1.0/(2.0*dx)

        ! Top and Bottom Boundary
        do j = 1, Ny
            dMx(1,j) = (M(2,j) - M(Nx,j))*dx_inv    ! Top
            dMx(Nx,j) = (M(1,j) - M(Nx-1,j))*dx_inv ! Bottom
        end do 

        ! Bulk (non-boundary nodes)
        do j = 1, Ny
            do i = 2, Nx-1
                dMx(i,j) = (M(i+1,j) - M(i-1,j))*dx_inv
            end do
        end do

    end subroutine dM_dx
    


    subroutine time_evolution(grid, grid_new, dQ, M, dt, Nx, Ny)
        ! Subroutine to perform the time evolution in accordance with
        ! the Cahn-Hilliard equation

        integer, intent(in) :: Nx, Ny
        real(real64), dimension(Nx, Ny), intent(in) :: grid, dQ
        real(real64), intent(in) :: M, dt
        real(real64), dimension(Nx, Ny), intent(out) :: grid_new
        integer :: i, j ! counters

        do i = 1, Nx
            do j = 1, Ny
                grid_new(i, j) = grid(i, j) + dt * M * dQ(i, j)
            end do
        end do

    end subroutine time_evolution

end module cahn_hilliard
