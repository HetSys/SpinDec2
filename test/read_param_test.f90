program main
    use input_params
    use iso_fortran_env

    implicit none

    integer  :: nx,ny,cint,random_seed,err
    character(len=128) :: cpi,cpo
    real(kind=REAL64) :: initial_conc, conc_std,m1,m2,k,bfe,t,delta_t,df_tol
    real(kind=real64),dimension(:),allocatable :: coeffs

    call read_params("../test/input_test.txt",initial_conc,conc_std,coeffs,nx,&
        ny,m1,m2,k,bfe,cint,cpi,cpo,t,delta_t,df_tol,random_seed,err)


    !print*,err
    if(err == -1) then
        print*, "There was an issue with the input file please check and try again"
        return
    end if

    print*,random_seed
end program main
