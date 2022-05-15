! Module containing subroutines for testing potentials
module test_potentials

    use iso_fortran_env
    use potentials

    implicit none

contains

    ! Subroutine that tests the bulk chemical potential
    ! @param c: concentration grid used for testing
    ! @param a: array storing coefficients
    ! @param expected: expected answer for bulk chemical potential
    subroutine test_bulk_potential(c, a, expected, test_num)

        real(real64), allocatable :: mu(:, :)
        real(real64), intent(in) :: c(:, :)
        real(real64), intent(in) :: a(:)
        real(real64), intent(in) :: expected(:, :)
        integer, intent(in) :: test_num
        logical :: res = .true.
        integer :: nx, ny, i, j

        nx = size(c, 1)
        ny = size(c, 2)

        allocate (mu(nx, ny))

        call bulk_potential(mu, c, a)

        ! Testing if resulting mu is equal to expected
        ! within certain tolerance
        do j = 1, ny
            do i = 1, nx
                if (abs(mu(i, j) - expected(i, j)) > 1e-3) then
                    res = .false.
                    exit
                end if
            end do
        end do

        if (res) then
            print *, 'Unit test #', test_num, ' for bulk potential succeeded.'
        else
            print *, 'Unit test #', test_num, ' for bulk potential failed.'
        end if

        deallocate (mu)

    end subroutine test_bulk_potential

    ! subroutine that tests the total chemical potential
    ! @param mu: bulk chemical potential grid used for testing
    ! @param c: concentration grid used for testing
    ! @param dx: spatial grid spacing in x direction used for testing
    ! @param dy: spatial grid spacing in y direction used for testing
    ! @param kappa: free energy gradient parameter used for testing
    ! @param expected: expected answer for total chemical potential
    subroutine test_total_potential(mu, c, dx, dy, Kappa, expected, test_num, analytical)

        real(real64), allocatable :: Q(:, :)
        real(real64), intent(in) :: c(:, :)
        real(real64), intent(in) :: mu(:, :)
        real(real64), intent(in) :: dx, dy, Kappa
        real(real64), intent(in) :: expected(:, :)
        integer, intent(in) :: test_num
        logical, intent(in) :: analytical
        logical :: res = .true.
        integer :: nx, ny, i, j,start_i,start_j,end_i,end_j

        nx = size(c, 1)
        ny = size(c, 2)

        allocate (Q(nx, ny))
        call total_potential(Q, mu, c, dx, dy, Kappa)

        if (analytical) then
            start_i = 2
            start_j = 2
            end_i = nx-1
            end_j = ny-1 
        else
            start_i = 1
            start_j = 1
            end_i = nx
            end_j = ny
        end if

        ! Testing if resulting Q is equal to expected
        ! within certain tolerance
        do j = start_j, end_j
            do i = start_i, end_i
                if (abs(Q(i, j) - expected(i, j)) > 1e-3) then
                    res = .false.
                    exit
                end if
            end do
        end do

        if (res) then
            print *, 'Unit test #', test_num, ' for total potential succeeded.'
        else
            print *, 'Unit test #', test_num, ' for total potential failed.'
        end if

        deallocate (Q)

    end subroutine test_total_potential

end module test_potentials
