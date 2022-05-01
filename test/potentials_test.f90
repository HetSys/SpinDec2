module test_potentials

    use iso_fortran_env
    use potentials
    
    implicit none

    contains
   
    subroutine test_bulk_potential(c,a,expected)
        
        real(real64),allocatable :: mu(:,:)
        real(real64),intent(in) :: c(:,:)
        real(real64),intent(in) :: a(:)
        real(real64),intent(in) :: expected(:,:)
        logical :: res = .true.
        integer :: nx,ny,i,j

        nx = size(c,1)
        ny = size(c,2)

        
        allocate(mu(nx,ny))

        call bulk_potential(mu,c,a)  

        do j=1,ny
            do i=1,nx
                if (abs(mu(i,j)-expected(i,j))>1e-3) then
                    res = .false.
                    exit
                end if        
            end do
        end do

        if (res) then
            print *, 'Unit test for bulk potential succeeded.'
        else
            print *, 'Unit test for bulk potential failed.'
        end if


        deallocate(mu)    
                
    end subroutine  

    subroutine test_total_potential(mu,c,dx,dy,kappa,expected)
        
        real(real64),allocatable :: Q(:,:)
        real(real64),intent(in)  :: c(:,:)
        real(real64),intent(in)  :: mu(:,:)
        real(real64),intent(in) :: dx,dy,kappa
        real(real64),intent(in) :: expected(:,:)
        logical :: res = .true.
        integer :: nx,ny,i,j

        nx = size(c,1)
        ny = size(c,2)

        allocate(Q(nx,ny))
        call total_potential(Q,mu,c,dx,dy,kappa)

        do j=1,ny
            do i=1,nx
                if (abs(Q(i,j)-expected(i,j))>1e-3) then
                    res = .false.
                    exit
                end if
            end do
        end do


        if (res) then
            print *, 'Unit test for total potential succeeded.'
        else
            print *, 'Unit test for total potential failed.'
            print *, Q
        end if

        deallocate(Q)

    end subroutine 
end module

program main

    use iso_fortran_env
    use test_potentials
    
    implicit none

    real(real64),dimension(3) :: a_1
    real(real64),dimension(3,3) :: c_1
    real(real64),dimension(3,3) :: expected_bulk_1,expected_total_1
    real(real64) :: dx,dy,kappa 

    print *, 'Started testing.'

    !Test 1 for bulk potential
    a_1 = (/1.0,2.0,3.0/)
    c_1 = reshape((/1.0,2.0,3.0,2.0,3.0,6.0,1.0,2.0,4.0/),shape(c_1))
    expected_bulk_1 = reshape((/8.0,14.0,20.0,14.0,20.0,38.0,8.0,14.0,26.0/),shape(c_1))
    call test_bulk_potential(c_1,a_1,expected_bulk_1)

    !Test 2 for bulk potential
    !call test_bulk_potential(c_2,a_2,expected_bulk_2)
    !Test 3 for bulk potential
    !call test_bulk_potential(c_3,a_3,expeceted_bulk_3)

    !Test 1 for total potential
    dx = 1.0
    dy = 1.0
    kappa = 2.0
    expected_total_1 = reshape((/0.0,12.0,18.0,8.0,20.0,62.0,-2.0,10.0,34.0/),shape(c_1))
    call test_total_potential(expected_bulk_1,c_1,dx,dx,kappa,expected_total_1)

    !Test 2 for total potential
    !call test_total_potential(mu_2,c_2,dx,dx,kappa,expected_total_2)
    !Test 3 for total potential
    !call test_total_potential(mu_3,c_3,dx,dx,kappa,expected_total_3)
 


end program



