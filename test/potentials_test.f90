module test_potentials

    use iso_fortran_env
    use potentials
    
    implicit none

    contains
   
    subroutine test_bulk_potential(c,a,expected)
        
        real(real64),allocatable :: mu(:,:)
        real(real64),intent(in) :: c(:,:)
        real(real64),intent(in) :: a(:)
        
        allocate(mu(size(c,1),size(c,2))
        call bulk_potential(mu,c,a)  
         
        res = abs(mu-expected)<1e-3

        if (res) then
            print *, 'Unit test for bulk potential succeeded.'
        else
            print *, 'Unit test for bulk potential failed.'
        end if


        deallocate(mu)    
                
    end subroutine  

    subroutine test_total_potential(mu,c,dx,dy,kappa,res,expected)
        
        real(real64),allocatable :: Q(:,:)
        real(real64),intent(in)  :: c(:,:)
        real(real64),intent(in) :: dx,dy,kappa

        allocate(Q(size(c,1),size(c,2))
        call total_potential(Q,mu,c,dx,dy,kappa)

        res = abs(Q-expected)<1e-3

        if (res) then
            print *, 'Unit test for bulk potential succeeded.'
        else
            print *, 'Unit test for bulk potential failed.'
        end if

        deallocate(Q)

    end subroutine 
end module

program main

    use iso_fortran_env
    use test_potentials
    
    implicit none

        

end program



