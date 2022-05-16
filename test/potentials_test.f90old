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

    end subroutine

    ! subroutine that tests the total chemical potential
    ! @param mu: bulk chemical potential grid used for testing
    ! @param c: concentration grid used for testing
    ! @param dx: spatial grid spacing in x direction used for testing
    ! @param dy: spatial grid spacing in y direction used for testing
    ! @param kappa: free energy gradient parameter used for testing
    ! @param expected: expected answer for total chemical potential
    subroutine test_total_potential(mu, c, dx, dy, Kappa, expected, test_num)

        real(real64), allocatable :: Q(:, :)
        real(real64), intent(in) :: c(:, :)
        real(real64), intent(in) :: mu(:, :)
        real(real64), intent(in) :: dx, dy, Kappa
        real(real64), intent(in) :: expected(:, :)
        integer, intent(in) :: test_num
        logical :: res = .true.
        integer :: nx, ny, i, j

        nx = size(c, 1)
        ny = size(c, 2)

        allocate (Q(nx, ny))
        call total_potential(Q, mu, c, dx, dy, Kappa)

        ! Testing if resulting Q is equal to expected
        ! within certain tolerance
        outer: do j = 1, ny
            do i = 1, nx
                if (abs(Q(i, j) - expected(i, j)) > 1e-3) then
                    res = .false.
                    print *, i,j
                    exit outer
                end if
            end do
        end do outer

        if (res) then
            print *, 'Unit test #', test_num, ' for total potential succeeded.'
        else
            print *, 'Unit test #', test_num, ' for total potential failed.'
        end if

        deallocate (Q)

    end subroutine
end module

! Program to run tests for the potentials
program main

    use iso_fortran_env
    use test_potentials

    implicit none

    real(real64), dimension(3) :: a_1
    real(real64), dimension(3, 3) :: c_1
    real(real64), dimension(3, 3) :: expected_bulk_1, expected_total_1
    real(real64) :: dx, dy, kappa
    real(real64), dimension(5) :: a_2
    real(real64), dimension(4, 4) :: c_2
    real(real64), dimension(4, 4) :: expected_bulk_2, expected_total_2
    integer :: test_num

    print *, 'Started testing.'

    !Test 1 for bulk potential
    a_1 = (/1.0, 2.0, 3.0/)
    c_1 = reshape((/1.0, 2.0, 3.0, 2.0, 3.0, 6.0, 1.0, 2.0, 4.0/), shape(c_1))
    expected_bulk_1 = reshape((/8.0, 14.0, 20.0, 14.0, 20.0, 38.0, 8.0, 14.0, 26.0/), shape(c_1))
    test_num = 1
    call test_bulk_potential(c_1, a_1, expected_bulk_1, test_num)

    !Test 2 for bulk potential
    a_2 = (/0.2, 0.5, 0.3, 1.4, 1.6/)
    c_2 = reshape((/0.16, 0.23, 0.32, 0.45, 0.35, 0.42, 0.71, 0.12, 0.24, 0.45, 0.64, 0.83, 0.22, 0.42, 0.18, 0.71/), shape(c_2))
    expected_bulk_2 = reshape((/0.7297, 0.9380, 1.3317, 2.2037, 1.4989, 1.9670, 5.3338, 0.6435, 0.9743, 2.2037, &
                                4.2820, 7.5508, 0.9034, 1.9670, 0.7814, 5.3338/), shape(c_2))
    test_num = 2
    call test_bulk_potential(c_2, a_2, expected_bulk_2, test_num)

    !Test 3 for bulk potential
    !call test_bulk_potential(c_3,a_3,expeceted_bulk_3)

    !Test 1 for total potential
    dx = 1.0
    dy = 1.0
    kappa = 2.0
    expected_total_1 = reshape((/0.0, 12.0, 18.0, 8.0, 20.0, 62.0, -2.0, 10.0, 34.0/), shape(c_1))
    test_num = 1
    call test_total_potential(expected_bulk_1, c_1, dx, dx, kappa, expected_total_1, test_num)

    !Test 2 for total potential
    expected_total_2 = reshape((/-0.4902, 0.1380, 0.7517, 3.1837, 2.4189, 1.8470, 8.0138, -3.0764, -0.8056, 2.3637, &
                                 5.0620, 10.7708, -0.3965, 3.1670, -1.9585, 7.6538/), shape(c_2))
    test_num = 2
    call test_total_potential(expected_bulk_2, c_2, dx, dx, kappa, expected_total_2, test_num)

    !Test 3 for total potential
    !call test_total_potential(mu_3,c_3,dx,dx,kappa,expected_total_3)

end program
