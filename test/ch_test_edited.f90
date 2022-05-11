module ch_test

    use iso_fortran_env
    use grid
    use cahn_hilliard
    use potentials

    implicit none

    contains

    subroutine test_del_Q(dQ,Q_test,dQ_expected,dx,dy,Nx,Ny)

        integer, intent(in) :: Nx, Ny
        real(real64), dimension(Nx,Ny), intent(in) :: Q_test, dQ_expected
        real(real64), dimension(Nx,Ny), intent(out) :: dQ
        real(real64), intent(in) :: dx, dy
        logical :: result = .TRUE.
        integer :: i,j ! counters

        call del_Q(dQ,Q_test,dx,dy,Nx,Ny)


        do i = 1, Nx
            do j = 1 , Ny
                ! print*, abs(dQ(i,j)-dQ_expected(i,j))
                if(abs(dQ(i,j)-dQ_expected(i,j)) > 1e-2) then
                    result = .FALSE.
                    exit ! do not continue checking
                end if
            end do
        end do

        if (result) then
            print*, "del_Q passes test"
        else
            print*, "del_Q does not pass test"
        end if

    end subroutine test_del_Q


    subroutine test_time_evolution(c_new,c_test,c_expected,Nx,Ny,dx,dy,dt,a,Kappa,M)

        integer, intent(in) :: Nx, Ny
        real(real64), dimension(Nx,Ny), intent(in) :: c_test, c_expected
        real(real64), intent(in) :: dx, dy, dt, Kappa, M
        real(real64), dimension(:), intent(in) :: a
        real(real64), dimension(Nx,Ny), intent(out) :: c_new
        real(real64), dimension(Nx,Ny) :: mu, Q, dQ
        logical :: result = .TRUE.
        integer :: i,j ! counters


        ! Get bulk chemical potential
        call bulk_potential(mu,c_test,a)
        ! Get total chemical potential
        call total_potential(Q,mu,c_test,dx,dy,Kappa)
        ! Get dQ
        call del_Q(dQ,Q,dx,dy,Nx,Ny)
        ! Get grid at next time step
        call time_evolution(c_test,c_new,dQ,M,dt,Nx,Ny)

        call del_Q(dQ,Q,dx,dy,Nx,Ny)

        do i = 1, Nx
            do j = 1 , Ny
                ! print *, abs(c_new(i,j)-c_expected(i,j))
                if(abs(c_new(i,j)-c_expected(i,j)) > 1e-3) then
                    result = .FALSE.
                    exit ! do not continue checking
                end if
            end do
        end do

        if (result) then
            print*, "time_evolution passes test"
        else
            print*, "time_evolution does not pass test"
        end if

    end subroutine test_time_evolution

end module ch_test
