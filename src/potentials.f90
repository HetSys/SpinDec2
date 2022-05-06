module potentials

    use iso_fortran_env

    implicit none

    contains      
    
    ! Subroutine to calculate bulk chemical potential
    ! @param mu: 2D grid to store bulk chemical potentials
    ! @param c: 2D concentration grid
    ! @param a: 1D array storing coefficients provided by user  
    subroutine bulk_potential(mu,c,a)

        real(real64) :: mu(:,:)
        real(real64), intent(in) :: c(:,:)
        real(real64), intent(in) :: a(:)
        integer :: nx,ny,n,i,j,k

        !Get mu grid size and number of user input coefficients
        nx = size(mu,1)
        ny = size(mu,2)
        n = size(a)

        !Initialize mu array with zero values
        mu = 0.0

        ! Loop filling the mu array
        do j=1,ny
            do i=1,nx
                do k=1,n-1
                    mu(i,j) = mu(i,j) + k*a(k+1)*c(i,j)**(k-1)
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
    subroutine total_potential(Q,mu,c,dx,dy,Kappa)

        real(real64), intent(in) :: mu(0:,0:)
        real(real64), intent(in) :: c(0:,0:)
        real(real64), intent(in) :: dx,dy,Kappa
        real(real64) :: Q(0:,0:)
        real(real64) :: x_lap,y_lap
        integer :: nx,ny,i,j

        ! Get Q grid size
        nx = size(Q,1)
        ny = size(Q,2)

        !Loop filling the Q array
        do j=0,ny-1
            do i=0,nx-1
                !X and Y double derivatives using finite differences
                !Peridic boundary conditions are taken care using mod function
                x_lap = (c(modulo(i+1,nx),j)-2*c(i,j)+c(modulo(i-1,nx),j))/(dx*dx)
                y_lap = (c(i,modulo(j+1,ny))-2*c(i,j)+c(i,modulo(j-1,ny)))/(dy*dy)
                Q(i,j) = mu(i,j) - Kappa*(x_lap + y_lap)
            end do
        end do

    end subroutine total_potential

end module potentials
