module spectral
    use fftw3

    use iso_fortran_env
    use potentials
    implicit none

contains

    subroutine spectral_method_iter(c_in, c_prev_in,a,dt,M,k, c_out,init)
        real(kind=real64) ,dimension(:,:), intent(in):: c_in,c_prev_in
        real(kind=real64) ,dimension(:,:), intent(inout):: c_out
        real(kind=real64) ,dimension(:,:), allocatable:: mu,k2,k4
        real(kind=real64),dimension(:), intent(in) :: a
        real(kind=real64) ,intent(in) :: dt,M,k
        integer :: i,j,init,dim12,dim22,stat
        complex(C_DOUBLE_COMPLEX),pointer,dimension(:,:) :: in,out_prev,in_prev, out,out_bulk,out_bulk_prev,ans
        integer, dimension(2) :: dims
        type(C_PTR) :: pin,pout,pout_prev,pout_bulk,pout_bulk_prev,plan,pans,plan_b
        real(kind=real64) :: norm, c_A
        real(kind=real64), parameter ::PI= 4*atan(1.0_real64)


        dims = shape(c_in)
        allocate(mu(dims(1),dims(2)))
        allocate(k2(dims(1),dims(2)))
        allocate(k4(dims(1),dims(2)))
        norm = sqrt(real(dims(1)*dims(2)))
        c_A = 0.1

        pin = fftw_alloc_complex(int(dims(1)*dims(2),C_SIZE_T))
        call c_f_pointer(pin, in, [dims(1),dims(2)])

        pout = fftw_alloc_complex(int(dims(1)*dims(2),C_SIZE_T))
        call c_f_pointer(pout, out, [dims(1),dims(2)])

        pout_prev = fftw_alloc_complex(int(dims(1)*dims(2),C_SIZE_T))
        call c_f_pointer(pout_prev, out_prev, [dims(1),dims(2)])

        pout_bulk = fftw_alloc_complex(int(dims(1)*dims(2),C_SIZE_T))
        call c_f_pointer(pout_bulk, out_bulk, [dims(1),dims(2)])

        pout_bulk_prev = fftw_alloc_complex(int(dims(1)*dims(2),C_SIZE_T))
        call c_f_pointer(pout_bulk_prev, out_bulk_prev, [dims(1),dims(2)])

        pans = fftw_alloc_complex(int(dims(1)*dims(2),C_SIZE_T))
        call c_f_pointer(pans, ans, [dims(1),dims(2)])


        in(:,:) = c_in(:,:)

        plan = fftw_plan_dft_2d(dims(1),dims(2),in,out,FFTW_FORWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(plan,in,out)
        !call fftw_destroy_plan(plan)

        if(init == 0) then
            in(:,:) = c_prev_in(:,:)

            !plan = fftw_plan_dft_2d(dims(1),dims(2),in,out_prev,FFTW_FORWARD,FFTW_ESTIMATE)
            call fftw_execute_dft(plan,in,out_prev)
            !call fftw_destroy_plan(plan)
        end if


        call bulk_potential(mu, c_in , a)
        in(:,:) = mu(:,:)

        !plan = fftw_plan_dft_2d(dims(1),dims(2),in,out_bulk,FFTW_FORWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(plan,in,out_bulk)
        !call fftw_destroy_plan(plan)


        if(init == 0) then
            call bulk_potential(mu, c_prev_in , a)
            in(:,:) = mu(:,:)

            !plan = fftw_plan_dft_2d(dims(1),dims(2),in,out_bulk_prev,FFTW_FORWARD,FFTW_ESTIMATE)
            call fftw_execute_dft(plan,in,out_bulk_prev)
            !call fftw_destroy_plan(plan)

        end if


        out = out
        out_bulk = out_bulk
        if(init == 0) then
            out_bulk_prev = out_bulk_prev
            out_prev = out_prev
        end if

        dim12 = int(real(dims(1))/2)
        dim22 = int(real(dims(2))/2)

        forall(i=1:dim12, j=1:dim22)
            k2(i,j)= 1+(-real(i-1)*real(i-1)-real(j-1)*real(j-1))*4*PI*PI/(dims(1)*dims(1))
        end forall

        forall(i= dim12+1:dims(1), j=dim22+1:dims(2))
            k2(i,j)= (-real(i-dims(1)-1)*real(i-dims(1)-1)-real(j-dims(2)-1)*real(j-dims(2)-1))*4*PI*PI/(dims(1)*dims(1))
        end forall

        forall(i= 1:dim12, j=dim22+1:dims(2))
            k2(i,j)= (-real(i-1)*real(i-1)-real(j-dims(2)-1)*real(j-dims(2)-1))*4*PI*PI/(dims(1)*dims(1))
        end forall

        forall(i= dim12+1:dims(1), j=1:dim22)
            k2(i,j)= (-real(i-dims(1)-1)*real(i-dims(1)-1)-real(j-1)*real(j-1))*4*PI*PI/(dims(1)*dims(1))
        end forall

        !open(unit=1, iostat=stat, file='ksq.dat', status='old')
        !if (stat == 0) close(1, status='delete')
        !
        !open(unit=1, iostat=stat, file='ksq.dat', status='new')
        !if (stat == 0) then
        !  do i=1,dims(1)
        !    do j=1,dims(2)
        !      write(1,'(F0.10)') k2(i,j)
        !    end do
        !  end do

        !  close(1)
        !end if


        !print*, k2


        forall(i=1:dims(1), j=1:dims(2))
            k4(i,j)=k2(i,j)*k2(i,j)
        end forall
        !print*, M
        !print*,k2
        if(init == 0) then
            ans = (2.0*dt/(3.0+2.0*dt*M*k*k4+2*dt*c_A*M*k2))*(2.0*M*k2*out_bulk-M*k2*&
                out_bulk_prev + M*k2*2*c_a*out - M*k2*c_a*out_prev + (4.0*out-out_prev)/(2.0*dt))
        else
            ans = dt*(-M*k*k4*out + M*k2*out_bulk)+out
        end if

        !print*,ans


        plan_b = fftw_plan_dft_2d(dims(1),dims(2),ans,in,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(plan_b,ans,in)
        call fftw_destroy_plan(plan_b)
        call fftw_destroy_plan(plan)

        c_out(:,:) = in(:,:)/(norm*norm)

        call fftw_free(pin)
        call fftw_free(pans)
        call fftw_free(pout)
        call fftw_free(pout_prev)
        call fftw_free(pout_bulk_prev)
        call fftw_free(pout_bulk)
        deallocate(k2)
        deallocate(k4)
        deallocate(mu)

    end subroutine spectral_method_iter



end module spectral
