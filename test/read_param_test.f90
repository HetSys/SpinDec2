program main
    use input_params
    use checkpointing
    use netcdf
    use iso_fortran_env

    implicit none

    integer  :: nx,ny,cint,random_seed,err,i,count,ncerr,use_input,current_time
    character(len=128) :: cpi,cpo
    real(kind=REAL64) :: initial_conc, conc_std,m1,m2,k,bfe,t,delta_t,df_tol
    real(kind=real64),dimension(:),allocatable :: coeffs
    real(kind=real64), dimension(:,:,:), allocatable :: datas
    real(kind=real64), dimension(:,:), allocatable :: mu
    real(kind=real64), dimension(:), allocatable :: ftot

    TYPE(run_data_type) :: run_data

    allocate(datas(1,1,1))
    allocate(mu(1,1))
    allocate(ftot(1))


    datas(1,1,1) = 1
    mu(1,1) = 1
    ftot(1) = 1
    current_time = 0
    count = 0

    !This is the way to call read_params. You have to input a file.
    !You must put all of the possible variables as an input. Optional ones will
    !default to -1 except the polynomial coeffs which will be an uninitialsed array
    !It will also require an intger that represents if an error has occured
    call read_params("../test/input_test.txt",initial_conc,conc_std,coeffs,nx,&
        ny,m1,m2,k,bfe,cint,cpi,cpo,t,delta_t,df_tol,random_seed,use_input,err)

    !If err is -1 then the input file was not read properly and the code needs to exit
    if(err == -1) then
        print*, "There was an issue with the input file please check and try again"
        stop
    end if

    !!This will be how the checkpoint input file will be read. IF it cannot be read the code will stop
    !!There is a flag use_input which will indicate whether to overwrite input or not
    if(cpi /= "") then
        call read_checkpoint_in(datas,mu,ftot, cpi,initial_conc,conc_std,coeffs,nx,&
            ny,m1,m2,k,bfe,cint,cpo,t,delta_t,df_tol,current_time,random_seed,use_input,ncerr)
        if(ncerr /= nf90_noerr) then
            print*, "There was an error reading the checkpoint file."
            stop
        end if
    end if




    !run_data%initial_conc = initial_conc
    !run_data%conc_std = conc_std
    !run_data%nx = nx
    !run_data%ny = ny
    !run_data%m1 = m1
    !run_data%k = k
    !run_data%bfe = bfe
    !run_data%cint = cint
    !run_data%t = t
    !run_data%time_step = delta_t
    !run_data%df_tol = df_tol
    !run_data%df_tol = df_tol
    !run_data%random_seed = random_seed
    !run_data%current_time =current_time

    !!This will be the format in which the checkpointing system will work

    do while(current_time < t/delta_t)
        if(count >= cint) then
            call write_checkpoint_file(datas,mu,ftot,coeffs,cpo,initial_conc,conc_std&
                ,nx,ny,m1,m2,k,bfe,Cint,t,delta_t,current_time,df_tol,random_seed,ncerr)
            if(ncerr /= nf90_noerr) then
                print*, "There was an error writing the checkpoint file."
                stop
            end if
            count = 0
        end if
        count = count +1
        current_time = current_time+1
        datas(1,1,1) = datas(1,1,1)+1
        print*, datas(1,1,1)

    end do




    print*,random_seed
end program main
