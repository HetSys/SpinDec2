module free_energy

    use iso_fortran_env

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

    subroutine bulk_free_energy(f_b,c,a)

        real(real64), intent(out), allocatable :: f_b(:,:)
        real(real64), intent(in) :: c(:,:)
        real(real64), intent(in) :: a(:)
        integer :: nx,ny,n
        integer :: i,j,k

        !Get f_b grid size and number of user input coefficients
        Nx = size(c,1)
        Ny = size(c,2)
        n = size(a)

        allocate(f_b(Nx,Ny))

        f_b = 0

        !Loop filling the f_b array
        do j=1,Ny
            do i=1,Nx
                do k=1,n
                    f_b(i,j) = f_b(i,j) + a(k)*c(i,j)**(k-1)
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

    subroutine total_free_energy(F, c, f_b, dx, dy, kappa)

        real(real64), intent(in) :: c(0:,0:)
        real(real64) :: f_b(0:,0:)
        integer :: Nx, Ny
        integer :: i,j
        real(real64), intent(out) :: F
        real(real64) :: dx, dy
        real(real64) :: kappa
        real(real64) :: P, xlap, ylap

        !Get mu grid size and number of user input coefficients
        Nx = size(c,1)
        Ny = size(c,2)

        F = 0.0

        !Loop filling the mu array
        do j=0,ny-1
             do i=0,nx-1

                xlap = (c(modulo(i+1,nx),j)-2*c(i,j)+c(modulo(i-1,nx),j))/(dx*dx)
                ylap = (c(i,modulo(j+1,ny))-2*c(i,j)+c(i,modulo(j-1,ny)))/(dy*dy)
                P = xlap + ylap

                F = F + (f_b(i,j) + 0.5*kappa*P*P)*dx*dy

                end do
            end do

    end subroutine



end module free_energy
