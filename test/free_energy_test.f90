! Module containing subroutines for testing free energies
module test_free_energy

    use iso_fortran_env
    use free_energy

    implicit none

    contains

    !******************************************************************************
    !> test_bulk_free
    !!
    !! subroutine to unit test bulk free energy
    !!
    !!@param c : concentration field
    !!@param a : polynomial coefficients
    !!@param expected : expected result for bulk free energy given a and c
    !!@param test_num : unit test number
    !******************************************************************************
    !
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

    !******************************************************************************
    !> test_total_free_energy
    !!
    !! subroutine to unit test total free energy
    !!
    !!@param c : concentration field
    !!@param f_b : bulk free energy field
    !!@param dx,dy : finite difference spacing
    !!@param kappa : interfacial term coefficient
    !!@param expected : expected result for bulk free energy given a and c
    !!@param test_num : unit test number
    !******************************************************************************
    subroutine test_total_free_energy(c, f_b, dx, dy, kappa,conc_halo, expected, test_num)


        real(real64), intent(in) :: c(:, :)
        real(real64), intent(in) :: f_b(:, :)
        real(real64), intent(in) :: conc_halo(:,:)
        real(real64), intent(in) :: dx, dy, kappa
        real(real64) :: expected
        integer, intent(in) :: test_num
        logical :: res = .true.
        integer :: nx, ny, i, j
        real(real64) :: F

        nx = size(c, 1)
        ny = size(c, 2)

        call total_free_energy(F, c, f_b, dx, dy, kappa,conc_halo)

        ! Testing if resulting F is equal to expected
        ! within certain tolerance
        if (abs(F - expected) > 1e-3) then
            res = .false.
        end if

        if (res) then
            print *, 'Unit test #', test_num, ' for total free energy succeeded.'
        else
            print *, 'Unit test #', test_num, ' for total free energy failed.'
        end if


    end subroutine test_total_free_energy

    !******************************************************************************
    !> free_energy_convergence
    !!
    !! subroutine to observe if F converges with reducing dx, dy
    !!
    !!@param c : concentration field
    !!@param f_b : bulk free energy field
    !!@param kappa : interfacial term coefficient
    !******************************************************************************


    subroutine free_energy_convergence(c, f_b, kappa,conc_halo)

        real(real64) :: dx, dy, delta
        real(real64), intent(in) :: kappa
        real(real64), intent(in) :: c(:, :)
        real(real64), intent(in) :: f_b(:, :)
        real(real64), intent(in) :: conc_halo(:,:)
        real(real64), dimension(0:10) :: F
        logical :: res = .false.
        integer :: i

        F = 0.0

        dx = 1
        dy = 1

        call total_free_energy(F(0), c, f_b, dx, dy, kappa,conc_halo)

        !observe whether F converges with increasing resoloution in x and y
        do i = 1, 10
            
            dx = dx/10
            dy = dy/10

            call total_free_energy(F(i), c, f_b, dx, dy, kappa,conc_halo)

            delta = abs(F(i) - F(i-1))

            if (delta < 1e-6) then
                res = .true.
                exit
             end if

        end do
        
        
        if (res) then
            print *, 'Convergence achieved'
        else
            print *, 'Convergence failed'
        end if

    end subroutine free_energy_convergence

end module test_free_energy
