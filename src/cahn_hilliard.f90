module cahn_hilliard

    use iso_fortran_env
    use grid
    use omp_lib

    implicit none

contains
 
    !> Subroutine to calculate diffusive mobility field.
    !!
    !!@param M ouput diffusive mobility field
    !!@param MA, MB user inputted atomic mobilities (one for each species)
    !!@param c0 user inputted mean of initial concentration field
    !!@param c concentration field at time step
    !!@param T temperature field
    !!@param problem character string that selects specific definition of M

    subroutine Mobility(M, MA, MB, EA, EB, c0, c, T, problem)

        real(real64), intent(inout):: M(:,:)
        real(real64), intent(in):: c(:,:), T(:,:)
        real(real64), intent(in):: MA, MB
        real(real64), intent(in):: EA, EB
        real(real64):: boltz
        real(real64):: c0
        real(real64), parameter:: R = 8.314
        integer:: Nx, Ny
        integer:: i, j
        character(len=*), intent(in):: problem  ! defines the set up of M


        Nx = size(c, 1)
        Ny = size(c, 2)

        select case (problem)

            !Temp-- Diffusive mobility dependent on c, and atomic mobilities
            !           are dependent on T
            case ("temp")

                !$omp parallel do default(shared) private(i, j, boltz)
                do j = 1, Ny
                    do i = 1, Nx
                        boltz = R*T(i, j)
                        M(i, j) = 1/boltz*exp(EA/boltz)*(1-c(i, j))
                        M(i, j) = M(i, j) + 1/boltz*exp(EB/boltz)*c(i, j)
                        M(i, j) = M(i, j)*c(i, j)*(1-c(i, j))
                    end do
                end do
                !$omp end parallel do

            !NonTemp-- Diffusive mobility still dependent on c, but atomic mobilities
            !           are independent on T
            case ("nontemp")

                !$omp parallel do default(shared) private(i, j)
                do j = 1, Ny
                    do i = 1, Nx
                        M(i, j) = (MA*(1-c(i, j)) + MB*c(i, j))*c(i, j)*(1-c(i, j))
                    end do
                end do
                !$omp end parallel do

            !Constant-Diffve Mobility is taken to be a mixing of the two constant
            !           mobilities, and is independent of c
            case ("constant")

                !$omp parallel do default(shared) private(i, j)
                do j = 1, Ny
                    do i = 1, Nx
                        M(i, j) = (MA*(1-c0) + MB*c0)*c0*(1-c0)
                    end do
                end do
                !$omp end parallel do

            case ("spectral")
                !$omp parallel do default(shared) private(i, j)
                do j = 1, Ny
                    do i = 1, Nx
                        M(i, j) = (MA*(1-c0) + MB*c0)*c0*(1-c0)
                    end do
                end do
                !$omp end parallel do
            ! condition to enforce one of the two setups
            case default
                print*, "PLEASE INPUT PROBLEM AS 'spectral', 'temp', 'nontemp' or 'constant'"
                stop "STOPPED"


            end select

    end subroutine


    !> Subroutine to calculate second dervative of Q
    !! using central differences
    !!
    !! @param dQ a 2D grid to store the second derivative of Q
    !! @param Q a 2D grid of total chemical potentials
    !! @param dx, dy the spatial grid spacings in x and y
    !! @param Nx, Ny grid size in x and y
    !! @param Q_halo array storing Q boundary neighbour values

    subroutine del_Q(dQ, Q, dx, dy, Nx, Ny, Q_halo)

        integer, intent(in):: Nx, Ny
        real(real64), dimension(Nx, Ny), intent(in):: Q
        real(real64), dimension(Nx, 4), intent(in):: Q_halo
        real(real64), dimension(Nx, Ny), intent(out):: dQ
        real(real64), intent(in):: dx, dy
        real(real64):: dx2, dy2
        real(real64):: der_x, der_y
        integer:: i, j  ! counters

        ! Pre-compute 1.0/(dx^2) and 1.0/(dy^2)
        dx2 = 1.0 / (dx*dx)
        dy2 = 1.0 / (dy*dy)

        ! Start computing derivatives
        ! Corner nodes
        ! Top-left
        der_x = (Q(2, 1) - 2.0*Q(1, 1) + Q_halo(1, left)) * dx2
        der_y = (Q(1, 2) - 2.0*Q(1, 1) + Q_halo(1, up)) * dy2
        dQ(1, 1) = (der_x+der_y)

        ! Bottom-left
        der_x = (Q(2, 1) - 2.0*Q(1, Ny) + Q_halo(1, left)) * dx2
        der_y = (Q_halo(1, down) - 2.0*Q(1, Ny) + Q(1, Ny-1)) * dy2
        dQ(1, Ny) = (der_x+der_y)

        ! Bottom-right
        der_x = (Q_halo(Ny, right) - 2.0*Q(Nx, Ny) + Q(Nx-1, Ny)) * dx2
        der_y = (Q_halo(Nx, down) - 2.0*Q(Nx, Ny) + Q(Nx, Ny-1)) * dy2
        dQ(Nx, Ny) = (der_x+der_y)

        ! Top-right
        der_x = (Q_halo(1, right) - 2.0*Q(Nx, 1) + Q(Nx-1, 1)) * dx2
        der_y = (Q(Nx, 2) - 2.0*Q(Nx, 1) + Q_halo(Nx, up)) * dy2
        dQ(Nx, 1) = (der_x+der_y)

        !print*, dQ(1, 1)

        ! Boundary nodes
        ! Top-j = 1
        !$omp parallel do default(shared) private(i, der_x, der_y)
        do i = 2, Nx-1
            der_y = (Q(i, 2) - 2.0*Q(i, 1) + Q_halo(i, up)) * dy2
            der_x = (Q(i+1, 1) - 2.0*Q(i, 1) + Q(i-1, 1)) * dx2
            dQ(i, 1) = (der_x+der_y)
            der_y = (Q_halo(i, down) - 2.0*Q(i, Ny) + Q(i, Ny-1)) * dy2
            der_x = (Q(i+1, Ny) - 2.0*Q(i, Ny) + Q(i-1, Ny)) * dx2
            dQ(i, Ny) = (der_x+der_y)
        end do
        !$omp end parallel do
        ! LHS-i = 1
        !$omp parallel do default(shared) private(j, der_x, der_y)
        do j = 2, Ny-1
            der_y = (Q(1, j+1) - 2.0*Q(1, j) + Q(1, j-1)) * dy2
            der_x = (Q(2, j) - 2.0*Q(1, j) + Q_halo(j, left)) * dx2
            dQ(1, j) = (der_x+der_y)
            der_y = (Q(Nx, j+1) - 2.0*Q(Nx, j) + Q(Nx, j-1)) * dy2
            der_x = (Q_halo(j, right) - 2.0*Q(Nx, j) + Q(Nx-1, j)) * dx2
            dQ(Nx, j) = (der_x+der_y)
        end do
        !$omp end parallel do
        ! Bulk (non-boundary) nodes
        !$omp parallel do default(shared) private(i, j, der_x, der_y)
        do i = 2, Nx-1
            do j = 2, Ny-1
                der_y = (Q(i, j+1) - 2*Q(i, j) + Q(i, j-1)) * dy2
                der_x = (Q(i+1, j) - 2*Q(i, j) + Q(i-1, j)) * dx2
                dQ(i, j) = (der_x+der_y)
            end do
        end do
        !$omp end parallel do

    end subroutine del_Q

    !> Subroutine to calculate the first dervative of Q
    !! with respect to x using central differences
    !!
    !! @param dQx a 2D grid to store the first derivative of Q
    !! @param Q a 2D grid of total chemical potentials
    !! @param dx, dy the spatial grid spacings in x and y
    !! @param Nx, Ny grid size in x and y
    !! @param Q_halo array storing Q boundary neighbour values

    subroutine dQ_dx(dQx, Q, dx, Nx, Ny, Q_halo)
        ! Subroutine to calculate dQ/dx
        ! Using central differences
        integer, intent(in):: Nx, Ny
        real(real64), dimension(Nx, Ny), intent(in):: Q
        real(real64), dimension(Nx, 4), intent(in):: Q_halo
        real(real64), intent(in):: dx
        real(real64), dimension(Nx, Ny), intent(out):: dQx
        real(real64):: dx_inv
        integer:: i, j  ! counters

        dx_inv = 1.0/(2.0*dx)

        ! LHS and RHS Boundary
        !$omp parallel do default(shared) private(j)
        do j = 1, Ny
            dQx(1, j) = (Q(2, j) - Q_halo(j, left))*dx_inv    ! Top
            dQx(Nx, j) = (Q_halo(j, right) - Q(Nx-1, j))*dx_inv  ! Bottom
            do i = 2, Nx-1
                dQx(i, j) = (Q(i+1, j) - Q(i-1, j))*dx_inv
            end do
        end do
        !$omp end parallel do

    end subroutine dQ_dx

    !> Subroutine to calculate the first dervative of Q
    !! with respect to y using central differences
    !!
    !! @param dQy a 2D grid to store the first derivative of Q
    !! @param Q a 2D grid of total chemical potentials
    !! @param dx, dy the spatial grid spacings in x and y
    !! @param Nx, Ny grid size in x and y
    !! @param Q_halo array storing Q boundary neighbour values

    subroutine dQ_dy(dQy, Q, dy, Nx, Ny, Q_halo)

        integer, intent(in):: Nx, Ny
        real(real64), dimension(Nx, Ny), intent(in):: Q
        real(real64), dimension(Nx, 4), intent(in):: Q_halo
        real(real64), intent(in):: dy
        real(real64), dimension(Nx, Ny), intent(out):: dQy
        real(real64):: dy_inv
        integer:: i, j  ! counters

        dy_inv = 1.0/(2.0*dy)

        ! Top and Bottom Boundary
        !$omp parallel do default(shared) private(i)
        do i = 1, Nx
            dQy(i, 1) = (Q(i, 2) - Q_halo(i, up))*dy_inv    ! LHS
            dQy(i, Ny) = (Q_halo(i, down) - Q(i, Ny-1))*dy_inv  ! RHS
        end do
        !$omp end parallel do
        ! Bulk (non-boundary nodes)
        !$omp parallel do default(shared) private(i, j)
        do j = 2, Ny-1
            do i = 1, Nx
                dQy(i, j) = (Q(i, j+1) - Q(i, j-1))*dy_inv
            end do
        end do
        !$omp end parallel do

    end subroutine dQ_dy


    !> Subroutine to calculate the first dervative of M
    !! with respect to y using central differences
    !!
    !! @param dMy a 2D grid to store the first derivative of M
    !! @param M a 2D grid of mobilities
    !! @param dx, dy the spatial grid spacings in x and y
    !! @param Nx, Ny grid size in x and y
    !! @param M_halo array storing M boundary neighbour values

    subroutine dM_dy(dMy, M, dy, Nx, Ny, M_halo)

        integer, intent(in):: Nx, Ny
        real(real64), dimension(Nx, Ny), intent(in):: M
        real(real64), dimension(Nx, 4), intent(in):: M_halo
        real(real64), intent(in):: dy
        real(real64), dimension(Nx, Ny), intent(out):: dMy
        real(real64):: dy_inv
        integer:: i, j  ! counters

        dy_inv = 1.0/(2.0*dy)

        ! Top and Bottom Boundary
        !$omp parallel do default(shared) private(i)
        do i = 1, Nx
            dMy(i, 1) = (M(i, 2) - M_halo(i, up))*dy_inv    ! Top
            dMy(i, Ny) = (M_halo(i, down) - M(i, Nx-1))*dy_inv  ! Bottom
        end do
        !$omp end parallel do
        ! Bulk (non-boundary nodes)
        !$omp parallel do default(shared) private(i, j)
        do j = 2, Ny-1
            do i = 1, Nx
                dMy(i, j) = (M(i, j+1) - M(i, j-1))*dy_inv
            end do
        end do
        !$omp end parallel do

    end subroutine dM_dy

    !> Subroutine to calculate the first dervative of Q
    !! with respect to y using central differences
    !!
    !! @param dMx a 2D grid to store the first derivative of Q
    !! @param M a 2D grid of mobilities
    !! @param dx, dy the spatial grid spacings in x and y
    !! @param Nx, Ny grid size in x and y
    !! @param M_halo array storing M boundary neighbour values

    subroutine dM_dx(dMx, M, dx, Nx, Ny, M_halo)
        ! Subroutine to calculate dM/dx
        ! Using central differences
        integer, intent(in):: Nx, Ny
        real(real64), dimension(Nx, Ny), intent(in):: M
        real(real64), dimension(Nx, 4), intent(in):: M_halo
        real(real64), intent(in):: dx
        real(real64), dimension(Nx, Ny), intent(out):: dMx
        real(real64):: dx_inv
        integer:: i, j  ! counters

        dx_inv = 1.0/(2.0*dx)
        ! LHS and RHS Boundary
        !$omp parallel do default(shared) private(j)
        do j = 1, Ny
            dMx(1, j) = (M(2, j) - M_halo(j, left))*dx_inv    ! LHS
            dMx(Nx, j) = (M_halo(j, right) - M(Nx-1, j))*dx_inv  ! RHS
            do i = 2, Nx-1
                dMx(i, j) = (M(i+1, j) - M(i-1, j))*dx_inv
            end do
        end do
        !$omp end parallel do
        !print*, dMx

    end subroutine dM_dx

    !> Subroutine to perform forward Euler (explicit) time integration
    !!
    !!@param c concentration field at previous step
    !!@param c_new concentration field at next step
    !!@param M diffusive mobility field
    !!@param dQ derivative of Q
    !!@param dx, dy grid spacings in x and y
    !!@param dt time-step
    !!@param Nx, Ny number of resoloutions in x and y

    subroutine time_evolution(grid, grid_new, dQ, M, dt, Nx, Ny)

        integer, intent(in):: Nx, Ny
        real(real64), dimension(Nx, Ny), intent(in):: grid, dQ
        real(real64), intent(in):: M, dt
        real(real64), dimension(Nx, Ny), intent(out):: grid_new
        integer:: i, j  ! counters

        do i = 1, Nx
            do j = 1, Ny
                grid_new(i, j) = grid(i, j) + dt*M * dQ(i, j)
            end do
        end do
    end subroutine time_evolution


    !> Subroutine to perform forward Euler (explicit) time integration
    !!
    !!@param c concentration field at previous step
    !!@param c_new concentration field at next step
    !!@param M diffusive mobility field
    !!@param Q total chemical potential
    !!@param dx, dy grid spacings in x and y
    !!@param dt time-step
    !!@param Nx, Ny number of resoloutions in x and y
    !!@param Q_halo array storing Q boundary neighbour values
    !!@param M_halo array storing M boundary neighbour values

    subroutine time_evolution_new(c, c_new, M, Q, dx, dy, dt, Nx, Ny, Q_halo, M_halo)
        integer, intent(in):: Nx, Ny
        real(real64), intent(in):: c(Nx, Ny), Q(Nx, Ny), M(Nx, Ny)
        real(real64), intent(in):: Q_halo(Nx, 4), M_halo(Nx, 4)
        real(real64):: dQ(Nx, Ny), dQx(Nx, Ny), dQy(Nx, Ny)
        real(real64):: dMx(Nx, Ny), dMy(Nx, Ny)
        real(real64), intent(out):: c_new(Nx, Ny)
        real(real64):: dx, dy, dt
        real(real64):: alpha, beta, xbeta, ybeta
        integer:: i, j  ! counters

        call del_Q(dQ, Q, dx, dy, Nx, Ny, Q_halo)
        call dQ_dx(dQx, Q, dx, Nx, Ny, Q_halo)
        call dM_dx(dMx, M, dx, Nx, Ny, M_halo)
        call dM_dy(dMy, M, dy, Nx, Ny, M_halo)
        call dQ_dy(dQy, Q, dy, Nx, Ny, Q_halo)

        !$omp parallel do default(shared) private(i, j, alpha, xbeta, ybeta, beta)
        do j = 1, Ny
            do i = 1, Nx
                alpha = M(i, j)*dQ(i, j)
                xbeta = dMx(i, j) * dQx(i, j)
                ybeta = dMy(i, j) * dQy(i, j)
                beta = xbeta+ybeta

                ! print*, alpha

                c_new(i, j) = c(i, j) + dt*(alpha+beta)
            end do
        end do
        !$omp end parallel do
    end subroutine time_evolution_new

end module cahn_hilliard
