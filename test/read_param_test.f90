program main
    use input_params
    use checkpointing
    use iso_fortran_env

    implicit none

    integer  :: nx,ny,cint,random_seed,err
    character(len=128) :: cpi,cpo
    real(kind=REAL64) :: initial_conc, conc_std,m1,m2,k,bfe,t,delta_t,df_tol
    real(kind=real64),dimension(:),allocatable :: coeffs
    real(kind=real64), dimension(:,:,:), allocatable :: data

    TYPE(run_data_type) :: run_data

    allocate(data(1,1,1))


    data(1,1,1) = 1

    !This is the way to call read_params. You have to input a file.
    !You must put all of the possible variables as an input. Optional ones will
    !default to -1 except the polynomial coeffs which will be an uninitialsed array
    !It will also require an intger that represents if an error has occured
    call read_params("../test/input_test.txt",initial_conc,conc_std,coeffs,nx,&
        ny,m1,m2,k,bfe,cint,cpi,cpo,t,delta_t,df_tol,random_seed,err)

    !If err is -1 then the input file was not read properly and the code needs to exit
    if(err == -1) then
        print*, "There was an issue with the input file please check and try again"
        return
    end if

    
    run_data%initial_conc = initial_conc
    run_data%conc_std = conc_std
    run_data%nx = nx
    run_data%ny = ny
    run_data%m1 = m1
    run_data%k = k
    run_data%bfe = bfe
    run_data%cint = cint
    run_data%t = t
    run_data%time_step = delta_t
    run_data%df_tol = df_tol
    run_data%df_tol = df_tol
    run_data%random_seed = random_seed
    run_data%current_time = 0

    call write_checkpoint_file(data,coeffs,cpo,run_data)


    print*,random_seed
end program main
