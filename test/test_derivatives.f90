module test_derivatives

    use iso_fortran_env
    use cahn_hilliard

    implicit none

    contains

    subroutine test_x_der(Nx, Ny)

        integer, intent(in) :: Nx, Ny
        integer ::  j, i
        real(real64) :: dx, dy
        real(real64) :: dx_inv, lh_bound, rh_bound
        real(kind=real64), dimension(Nx, Ny) :: test_in,test_out
        real(kind=real64), dimension(Nx, Ny) :: dQx, dMx, dQy, dMy
        real(kind=real64), dimension(Nx) :: x
        logical :: res = .true.
        real(real64) :: L
        REAL(REAL64) :: pi = 3.14159265_REAL64


        x = 0

        test_in = 0

        test_out = 0

        L = 1.0
        dx = L/(Nx-1)

        dx_inv = 1/(2.0*dx)
     
        do i = 1, Nx
            do j = 1, Ny

                x(i) = (i-1)*dx
                test_in(i,j) = sin(2*pi*x(i))
                test_out(i,j) = 2.0*pi*cos(2*pi*x(i))

            end do
        end do
    
        call dQ_dx(dQx,test_in,dx,Nx,Ny)

        ! test internal, should be the same as analytic

        ! test boundary works, forced to match c(0,j) to c(Nx,j)
        ! and c(Nx+1, j) to c(1, j)
        
        ! Testing if resulting dQx is equal to expected
        ! within certain tolerance
        do j = 1, Ny
            do i = 2, Nx-1
                if (abs(dQx(i, j) - test_out(i, j)) > 5e-3) then
                    res = .false.
                    exit
                end if
            end do
        end do

        if (res) then
            print *, 'Analytic test for dQ_dx succeeded.'
        else
            print *, 'Analytic test for dQ_dx failed.'
        end if

        res = .true.

        do j = 1, Ny
            lh_bound = (test_in(2,j) - test_in(Nx, j))*dx_inv
            rh_bound = (test_in(1,j) - test_in(Nx-1, j))*dx_inv

            if (abs(lh_bound - dQx(1,j)) > 1e-6) then
                res = .false.
                exit
            end if

            if (abs(rh_bound - dQx(Nx,j)) > 1e-6) then
                res = .false.
                exit
            end if

            if (abs(rh_bound - lh_bound) > 1e-6) then
                res = .false.
                exit
            end if

        end do

        if (res) then
            print *, 'dQ_dx PBC test succeded'
        else
            print *, 'dQ_dx PBC test Failed'
        end if

        ! dM/dx tests

        call dM_dx(dMx,test_in,dx,Nx,Ny)

        ! test internal, should be the same as analytic

        ! test boundary works, forced to match c(0,j) to c(Nx,j)
        ! and c(Nx+1, j) to c(1, j)
        
        ! Testing if resulting dMx is equal to expected
        ! within certain tolerance
        do j = 1, Ny
            do i = 2, Nx-1
                if (abs(dMx(i, j) - test_out(i, j)) > 5e-3) then
                    res = .false.
                    exit
                end if
            end do
        end do

        if (res) then
            print *, 'Analytic test for dM_dx succeeded.'
        else
            print *, 'Analytic test for dM_dx failed.'
        end if

        res = .true.

        do j = 1, Ny
            lh_bound = (test_in(2,j) - test_in(Nx, j))*dx_inv
            rh_bound = (test_in(1,j) - test_in(Nx-1, j))*dx_inv

            if (abs(lh_bound - dMx(1,j)) > 1e-6) then
                res = .false.
                exit
            end if

            if (abs(rh_bound - dMx(Nx,j)) > 1e-6) then
                res = .false.
                exit
            end if

            if (abs(rh_bound - lh_bound) > 1e-6) then
                res = .false.
                exit
            end if

        end do

        if (res) then
            print *, 'dM_dx PBC test succeded'
        else
            print *, 'dM_dx PBC test Failed'
        end if

        !test y derivatives if the grid give zero everywhere
        dy = dx

        call dQ_dy(dQy,test_in,dy,Nx,Ny)
        call dM_dy(dMy,test_in,dy,Nx,Ny)

        res = .true.

        do j = 1, Ny
            do i = 1, Nx
                if (abs(dMy(i, j)) > 1e-12) then
                    res = .false.
                    exit
                end if
                if (abs(dQy(i, j)) > 1e-12) then
                    res = .false.
                    exit
                end if
            end do
        end do

        if (res) then
            print *, 'first x deriavtive subroutines passed all tests'
        else
            print *, 'ERROR: y derivatve occurs in x derivative only subroutine'
        end if


    end subroutine test_x_der

    subroutine test_y_der(Nx, Ny)

        integer, intent(in) :: Nx, Ny
        integer ::  j, i
        real(real64) :: dy, dx
        real(real64) :: dy_inv, lh_bound, rh_bound
        real(kind=real64), dimension(Nx, Ny) :: test_in,test_out
        real(kind=real64), dimension(Nx, Ny) :: dQy, dMy, dQx, dMx
        real(kind=real64), dimension(Ny) :: y
        logical :: res = .true.
        real(real64) :: L
        REAL(REAL64) :: pi = 3.14159265_REAL64


        y = 0

        test_in = 0

        test_out = 0

        L = 1.0
        dy = L/(Nx-1)

        dy_inv = 1/(2.0*dy)
     
        do j = 1, Ny
            do i = 1, Nx

                y(j) = (j-1)*dy
                test_in(i,j) = sin(2*pi*y(j))
                test_out(i,j) = 2.0*pi*cos(2*pi*y(j))

            end do
        end do
    
        call dQ_dy(dQy,test_in,dy,Nx,Ny)

        ! test internal, should be the same as analytic

        ! test boundary works, forced to match c(0,j) to c(Nx,j)
        ! and c(Nx+1, j) to c(1, j)
        
        ! Testing if resulting dQy is equal to expected
        ! within certain tolerance
        do j = 2, Ny-1
            do i = 1, Nx
                if (abs(dQy(i, j) - test_out(i, j)) > 5e-3) then
                    res = .false.
                    exit
                end if
            end do
        end do

        if (res) then
            print *, 'Analytic test for dQ_dy succeeded.'
        else
            print *, 'Analytic test for dQ_dy failed.'
        end if

        res = .true.

        do i = 1, Nx
            lh_bound = (test_in(i,2) - test_in(i, Ny))*dy_inv
            rh_bound = (test_in(i,1) - test_in(i, Ny-1))*dy_inv

            if (abs(lh_bound - dQy(i,1)) > 1e-6) then
                res = .false.
                exit
            end if

            if (abs(rh_bound - dQy(i,Ny)) > 1e-6) then
                res = .false.
                exit
            end if

            if (abs(rh_bound - lh_bound) > 1e-6) then
                res = .false.
                exit
            end if

        end do

        if (res) then
            print *, 'dQ_dy PBC test succeded'
        else
            print *, 'dQ_dy PBC test Failed'
        end if

        ! dM/dx tests

        call dM_dy(dMy,test_in,dy,Nx,Ny)

        ! test internal, should be the same as analytic

        ! test boundary works, forced to match c(0,j) to c(Nx,j)
        ! and c(Nx+1, j) to c(1, j)
        
        ! Testing if resulting dMx is equal to expected
        ! within certain tolerance
        do j = 2, Ny-1
            do i = 1, Nx
                if (abs(dMy(i, j) - test_out(i, j)) > 5e-3) then
                    res = .false.
                    exit
                end if
            end do
        end do

        if (res) then
            print *, 'Analytic test for dM_dy succeeded.'
        else
            print *, 'Analytic test for dM_dy failed.'
        end if

        res = .true.

        do i = 1, Nx
            lh_bound = (test_in(i,2) - test_in(i, Ny))*dy_inv
            rh_bound = (test_in(i,1) - test_in(i, Ny-1))*dy_inv

            if (abs(lh_bound - dMy(i,1)) > 1e-6) then
                res = .false.
                exit
            end if

            if (abs(rh_bound - dMy(i,Ny)) > 1e-6) then
                res = .false.
                exit
            end if

            if (abs(rh_bound - lh_bound) > 1e-6) then
                res = .false.
                exit
            end if

        end do

        if (res) then
            print *, 'dM_dy PBC test succeded'
        else
            print *, 'dM_dy PBC test Failed'
        end if

        !test y derivatives if the grid give zero everywhere
        dx = dy

        call dQ_dx(dQx,test_in,dx,Nx,Ny)
        call dM_dx(dMx,test_in,dx,Nx,Ny)

        res = .true.

        do j = 1, Ny
            do i = 1, Nx
                if (abs(dMx(i, j)) > 1e-12) then
                    res = .false.
                    exit
                end if
                if (abs(dQx(i, j)) > 1e-12) then
                    res = .false.
                    exit
                end if
            end do
        end do

        if (res) then
            print *, 'first y deriavtive subroutines passed all tests'
        else
            print *, 'ERROR: x derivatve occurs in y derivative only subroutine'
        end if


    end subroutine test_y_der


end module test_derivatives


