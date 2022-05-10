! Module containing subroutines for testing free energies
module test_free_energy

    use iso_fortran_env
    use free_energy

    implicit none

    contains

    ! Subroutine that tests the bulk free energy
    ! @param c: concentration grid used for testing
    ! @param a: array storing coefficients
    ! @param expected: expected answer for bulk chemical potential
    ! @param test_num : test number
    subroutine test_bulk_free(c, a, expected, test_num)

        real(real64), allocatable :: f_b(:, :)
        real(real64), intent(in) :: c(:, :)
        real(real64), intent(in) :: a(:)
        real(real64), intent(in) :: expected(:, :)
        integer, intent(in) :: test_num
        logical :: res = .true.
        integer :: nx, ny, i, j

        nx = size(c, 1)
        ny = size(c, 2)

        allocate (f_b(nx, ny))

        call bulk_free_energy(f_b, c, a)

        ! Testing if resulting f_b is equal to expected
        ! within certain tolerance
        do j = 1, ny
            do i = 1, nx
                if (abs(f_b(i, j) - expected(i, j)) > 1e-3) then
                    res = .false.
                    exit
                end if
            end do
        end do

        if (res) then
            print *, 'Unit test #', test_num, ' for bulk free energy succeeded.'
        else
            print *, 'Unit test #', test_num, ' for free energy failed.'
        end if

        deallocate (f_b)

    end subroutine test_bulk_free

    ! subroutine that tests the total chemical potential
    ! @param f_b: bulk free energy grid used for testing
    ! @param c: concentration grid used for testing
    ! @param dx: spatial grid spacing in x direction used for testing
    ! @param dy: spatial grid spacing in y direction used for testing
    ! @param kappa: free energy gradient parameter used for testing
    ! @param expected: expected answer for total chemical potential
    ! @param test_num : test number
    subroutine test_total_free_energy(c, f_b, dx, dy, kappa, expected, test_num)


        real(real64), intent(in) :: c(:, :)
        real(real64), intent(in) :: f_b(:, :)
        real(real64), intent(in) :: dx, dy, kappa
        real(real64) :: expected
        integer, intent(in) :: test_num
        logical :: res = .true.
        integer :: nx, ny, i, j
        real(real64) :: F

        nx = size(c, 1)
        ny = size(c, 2)

        call total_free_energy(F, c, f_b, dx, dy, kappa)

        ! Testing if resulting F is equal to expected
        ! within certain tolerance
        do j = 1, ny
            do i = 1, nx
                if (abs(F - expected) > 1e-3) then
                    res = .false.
                    exit
                end if
            end do
        end do

        if (res) then
            print *, 'Unit test #', test_num, ' for total free energy succeeded.'
        else
            print *, 'Unit test #', test_num, ' for total free energy failed.'
        end if


    end subroutine test_total_free_energy
end module test_free_energy