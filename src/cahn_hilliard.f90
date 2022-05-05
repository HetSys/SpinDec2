module cahn_hilliard

    use iso_fortran_env

    implicit none

    contains

    subroutine del_Q(dQ,Q,M,dx,dy,Nx,Ny)
        ! Subroutine to perform spatial dervative of Q
        ! dQ is a 2D grid to store the derivatives of Q
        ! Q is the 2D grid of total chemical potentials
        ! M is the mobility of the species
        ! dx and dy are the spatial grid spacings in x and y

        implicit none

        integer, intent(in) :: Nx, Ny
        real(real64), dimension(Nx,Ny), intent(in) :: Q
        real(real64), dimension(Nx,Ny), intent(out) :: dQ
        real(real64), intent(in) :: M, dx, dy
        real(real64) :: dx2, dy2
        real(real64) :: der_x, der_y
        integer :: i,j ! counters

        ! Pre-compute 1.0/(dx^2) and 1.0/(dy^2)
        dx2 = 1.0/(dx*dx)
        dy2 = 1.0/(dy*dy)

        ! Start computing derivatives
        ! Corner nodes
        ! Top-left
        der_x = (Q(2,1)-2.0*Q(1,1)+Q(Nx,1))*dx2
        der_y = (Q(1,2)-2.0*Q(1,1)+Q(1,Ny))*dy2
        dQ(1,1) = M*(der_x + der_y)

        ! Bottom-left
        der_x = (Q(1,1)-2.0*Q(Nx,1)+Q(Nx-1,1))*dx2
        der_y = (Q(Nx,2)-2.0*Q(Nx,1)+Q(Nx,Ny))*dy2
        dQ(Nx,1) = M*(der_x + der_y)

        ! Bottom-right
        der_x = (Q(1,Ny)-2.0*Q(Nx,Ny)+Q(Nx-1,Ny))*dx2
        der_y = (Q(Nx,1)-2.0*Q(Nx,Ny)+Q(Nx,Ny-1))*dy2
        dQ(Nx,Ny) = M*(der_x + der_y)

        ! Top-right
        der_x = (Q(2,Ny)-2.0*Q(1,Ny)+Q(Nx,Ny))*dx2
        der_y = (Q(1,1)-2.0*Q(1,Ny)+Q(1,Ny-1))*dy2
        dQ(1,Ny) = M*(der_x + der_y)

        ! Boundary nodes
        ! LHS - j = 1
        do i = 2, Nx-1
            der_x = (Q(i+1,1)-2.0*Q(i,1)+Q(i-1,1))*dx2
            der_y = (Q(i,2)-2.0*Q(i,1)+Q(i,Ny))*dy2
            dQ(i,1) = M*(der_x + der_y)
        end do

        ! RHS - j = Ny
        do i = 2, Nx-1
            der_x = (Q(i+1,Ny)-2.0*Q(i,Ny)+Q(i-1,Ny))*dx2
            der_y = (Q(i,1)-2.0*Q(i,Ny)+Q(i,Ny-1))*dy2
            dQ(i,Ny) = M*(der_x + der_y)
        end do

        ! Top - i = 1
        do j = 2, Ny-1
            der_x = (Q(2,j)-2.0*Q(1,j)+Q(Nx,j))*dx2 
            der_y = (Q(1,j+1)-2.0*Q(1,j)+Q(1,j-1))*dy2
            dQ(1,j) = M*(der_x + der_y)
        end do

        ! Bottom - i = Nx
        do j = 2, Ny-1
            der_x = (Q(1,j)-2.0*Q(Nx,j)+Q(Nx-1,j))*dx2 
            der_y = (Q(Nx,j+1)-2.0*Q(Nx,j)+Q(Nx,j-1))*dy2
            dQ(Nx,j) = M*(der_x + der_y)
        end do

        ! Bulk (non-boundary) nodes
        do i = 2, Nx-1
            do j = 2, Ny-1
                der_x = (Q(i+1,j)-2*Q(i,j)+Q(i-1,j))*dx2
                der_y = (Q(i,j+1)-2*Q(i,j)+Q(i,j-1))*dy2
                dQ(i,j) = M*(der_x + der_y)
            end do
        end do
    
    end subroutine del_Q

    subroutine time_evolution(grid,grid_new,dQ,dt,Nx,Ny)
        ! Subroutine to perform the time evolution in accordance with
        ! the Cahn-Hilliard equation

        implicit none

        integer, intent(in) :: Nx, Ny
        real(real64), dimension(Nx,Ny), intent(in) :: grid, dQ
        real(real64), intent(in) :: dt
        real(real64), dimension(Nx,Ny), intent(out) :: grid_new
        integer :: i,j ! counters

        do i = 1, Nx
            do j = 1, Ny
                grid_new(i,j) = grid(i,j) + dt*dQ(i,j)
            end do
        end do

    end subroutine time_evolution

end module cahn_hilliard