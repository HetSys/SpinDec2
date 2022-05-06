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
    subroutine test_bulk_potential(c,a,expected,test_num)

        real(real64),allocatable :: mu(:,:)
        real(real64),intent(in) :: c(:,:)
        real(real64),intent(in) :: a(:)
        real(real64),intent(in) :: expected(:,:)
        integer,intent(in) :: test_num
        logical :: res = .true.
        integer :: nx,ny,i,j

        nx = size(c,1)
        ny = size(c,2)


        allocate(mu(nx,ny))

        call bulk_potential(mu,c,a)

        ! Testing if resulting mu is equal to expected
        ! within certain tolerance
        do j=1,ny
            do i=1,nx
                if (abs(mu(i,j)-expected(i,j))>1e-3) then
                    res = .false.
                    exit
                end if
            end do
        end do

        if (res) then
            print *, 'Unit test #',test_num,' for bulk potential succeeded.'
        else
            print *, 'Unit test #',test_num,' for bulk potential failed.'
        end if


        deallocate(mu)

    end subroutine test_bulk_potential

    ! subroutine that tests the total chemical potential
    ! @param mu: bulk chemical potential grid used for testing
    ! @param c: concentration grid used for testing
    ! @param dx: spatial grid spacing in x direction used for testing
    ! @param dy: spatial grid spacing in y direction used for testing
    ! @param kappa: free energy gradient parameter used for testing
    ! @param expected: expected answer for total chemical potential
    subroutine test_total_potential(mu,c,dx,dy,Kappa,expected,test_num)

        real(real64),allocatable :: Q(:,:)
        real(real64),intent(in)  :: c(:,:)
        real(real64),intent(in)  :: mu(:,:)
        real(real64),intent(in) :: dx,dy,Kappa
        real(real64),intent(in) :: expected(:,:)
        integer,intent(in) :: test_num
        logical :: res = .true.
        integer :: nx,ny,i,j

        nx = size(c,1)
        ny = size(c,2)

        allocate(Q(nx,ny))
        call total_potential(Q,mu,c,dx,dy,Kappa)

        ! Testing if resulting Q is equal to expected
        ! within certain tolerance
        do j=1,ny
            do i=1,nx
                if (abs(Q(i,j)-expected(i,j))>1e-3) then
                    res = .false.
                    exit
                end if
            end do
        end do


        if (res) then
            print *, 'Unit test #',test_num,' for total potential succeeded.'
        else
            print *, 'Unit test #',test_num,' for total potential failed.'
        end if

        deallocate(Q)

    end subroutine test_total_potential
end module test_potentials
