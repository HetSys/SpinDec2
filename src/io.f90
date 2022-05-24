module io

    use iso_fortran_env
    use netcdf

    implicit none

contains

    !************************************************************************
    !> check
    !!
    !! subroutine to perform error status check
    !!
    !!@param status : status code for error
    !************************************************************************
    !
    subroutine check(status)

        integer, intent(in) :: status

        if (status /= nf90_noerr) then
            print *, trim(nf90_strerror(status))
            stop ""
        end if

    end subroutine check

    !******************************************************************************
    !> write_netcdf
    !!
    !! subroutine to perform writing ouput file (netcdf)
    !!
    !!@param filename : filename to store data in
    !!@param c : order parameter field over all time
    !!@param F_tot : total free energy over time
    !!@param a : polynomial coefficients for f(c)
    !!@param Nc : number of coefficients
    !!@params Nx, Ny, Nt : number of  discretizations  in x, y and time
    !!@param dt : forward Euler time-step
    !!@param c0 : Initialaverage concentration
    !!@params MA, MB : Atomic mobilities of species A and B
    !!@param kappa : gradient term coefficient
    !*****************************************************************************
    !
    subroutine write_netcdf(c, F_tot, a, Nc, Nx, Ny, Nt, dt, c0, MA, MB, kappa)

        integer, intent(in) :: Nx, Ny, Nt, Nc
        real(kind=real64), dimension(Nx, Ny, Nt), intent(in) :: c
        real(kind=real64), dimension(Nt), intent(in) :: F_tot
        real(kind=real64), dimension(Nc), intent(in) :: a
        integer, dimension(3) :: size_c
        integer, dimension(3) :: c_dim_ids
        integer :: size_F_tot, size_a
        integer :: F_tot_dim_ids, a_dim_ids
        integer :: c_var_id, F_tot_var_id, a_var_id
        character(len=*), dimension(3), parameter :: c_dims = (/"c_x", "c_y", "c_t"/)
        character(len=*), parameter :: F_tot_dims = "F_t", a_dims = "a_i"
        character(len=*), parameter :: filename = 'CH_output2.nc'
        real(kind=real64), intent(in) :: dt, MA, MB, kappa, c0
        integer :: k, file_id

        !define size of arrays
        size_c = shape(c)
        size_F_tot = size(F_tot)
        size_a = size(a)

        !Create the file, overwriting if it exists
        call check(nf90_create(filename, NF90_CLOBBER, file_id))

        !write in the dimensions for the variables
        do k = 1, 3
            call check(nf90_def_dim(file_id, c_dims(k), size_c(k), c_dim_ids(k)))
        end do

        call check(nf90_def_dim(file_id, F_tot_dims, size_F_tot, F_tot_dim_ids))
        call check(nf90_def_dim(file_id, a_dims, size_a, a_dim_ids))
        print *, "Dimensions Written"

        !Define variable type, matching the arrays
        call check(nf90_def_var(file_id, "c", NF90_DOUBLE, c_dim_ids, c_var_id))
        call check(nf90_def_var(file_id, "F_tot", NF90_DOUBLE, F_tot_dim_ids, F_tot_var_id))
        call check(nf90_def_var(file_id, "coeffs", NF90_DOUBLE, a_dim_ids, a_var_id))
        print *, "Variable types defined"

        ! Setting up meta variables of simulation

        !Time step
        call check(nf90_put_att(file_id, NF90_GLOBAL, "dt", dt))

        !Number of time iterations
        call check(nf90_put_att(file_id, NF90_GLOBAL, "Nt", Nt))

        !Grid size (x)
        call check(nf90_put_att(file_id, NF90_GLOBAL, "Nx", Nx))

        !Grid size (y)
        call check(nf90_put_att(file_id, NF90_GLOBAL, "Ny", Ny))

        !Mobility
        call check(nf90_put_att(file_id, NF90_GLOBAL, "MA", MA))
        call check(nf90_put_att(file_id, NF90_GLOBAL, "MB", MB))

        !gradient term coefficient
        call check(nf90_put_att(file_id, NF90_GLOBAL, "kappa", kappa))

        !initial average concentration
        call check(nf90_put_att(file_id, NF90_GLOBAL, "c0", c0))

        !Finish defining metadata
        call check(nf90_enddef(file_id))

        ! Actually write the variables
        call check(nf90_put_var(file_id, c_var_id, c))
        call check(nf90_put_var(file_id, F_tot_var_id, F_tot))
        call check(nf90_put_var(file_id, a_var_id, a))

        !close the file
        call check(nf90_close(file_id))

        print *, "Sucessfully written netcdf file"

    end subroutine write_netcdf


    subroutine write_netcdf_parts_setup(c, F_tot, a, Nc, Nx, Ny, Nt, dt, c0, MA, MB, kappa,file_id)

        integer, intent(in) :: Nx, Ny, Nt, Nc
        real(kind=real64), dimension(Nx, Ny, Nt), intent(in) :: c
        real(kind=real64), dimension(Nt), intent(in) :: F_tot
        real(kind=real64), dimension(Nc), intent(in) :: a
        integer(4), dimension(3) :: size_c
        integer, dimension(3) :: c_dim_ids
        integer :: size_F_tot, size_a
        integer :: F_tot_dim_ids, a_dim_ids
        integer :: c_var_id, F_tot_var_id, a_var_id
        character(len=*), dimension(3), parameter :: c_dims = (/"c_x", "c_y", "c_t"/)
        character(len=*), parameter :: F_tot_dims = "F_t", a_dims = "a_i"
        character(len=*), parameter :: filename = 'CH_output.nc'
        real(kind=real64), intent(in) :: dt, MA, MB, kappa, c0
        integer :: k
        integer, intent(out) :: file_id

        !define size of arrays
        size_c = shape(c)
        size_c(3) = NF90_UNLIMITED
        !print*, NF90_UNLIMITED
        size_F_tot = size(F_tot)
        size_a = size(a)

        !Create the file, overwriting if it exists
        call check(nf90_create(filename, NF90_CLOBBER, file_id))

        !write in the dimensions for the variables
        do k = 1, 3
            call check(nf90_def_dim(file_id, c_dims(k), size_c(k), c_dim_ids(k)))
        end do

        call check(nf90_def_dim(file_id, F_tot_dims, size_F_tot, F_tot_dim_ids))
        call check(nf90_def_dim(file_id, a_dims, size_a, a_dim_ids))
        print *, "Dimensions Written"

        !Define variable type, matching the arrays
        call check(nf90_def_var(file_id, "c", NF90_DOUBLE, c_dim_ids, c_var_id))
        call check(nf90_def_var(file_id, "F_tot", NF90_DOUBLE, F_tot_dim_ids, F_tot_var_id))
        call check(nf90_def_var(file_id, "coeffs", NF90_DOUBLE, a_dim_ids, a_var_id))
        print *, "Variable types defined"

        ! Setting up meta variables of simulation

        !Time step
        call check(nf90_put_att(file_id, NF90_GLOBAL, "dt", dt))

        !Number of time iterations
        call check(nf90_put_att(file_id, NF90_GLOBAL, "Nt", Nt))

        !Grid size (x)
        call check(nf90_put_att(file_id, NF90_GLOBAL, "Nx", Nx))

        !Grid size (y)
        call check(nf90_put_att(file_id, NF90_GLOBAL, "Ny", Ny))

        !Mobility
        call check(nf90_put_att(file_id, NF90_GLOBAL, "MA", MA))
        call check(nf90_put_att(file_id, NF90_GLOBAL, "MB", MB))

        !gradient term coefficient
        call check(nf90_put_att(file_id, NF90_GLOBAL, "kappa", kappa))

        !initial average concentration
        call check(nf90_put_att(file_id, NF90_GLOBAL, "c0", c0))

        !Finish defining metadata
        call check(nf90_enddef(file_id))

        call check(nf90_put_var(file_id, a_var_id, a))

        !close the file
        !call check(nf90_close(file_id))

        print *, "Sucessfully setup netcdf file"

    end subroutine write_netcdf_parts_setup


    subroutine write_netcdf_parts(c, F_tot,PT,file_id)
        real(kind=real64), dimension(:, :,:), intent(in) :: c
        real(kind=real64),dimension(:), intent(in) :: F_tot
        character(len=*), parameter :: filename = 'CH_output.nc'
        integer :: k ,ierr
        integer, dimension(2) :: var_ids
        integer,intent(in) :: PT
        integer, intent(inout) :: file_id


        ierr = nf90_inq_varid(file_id, "F_tot", var_ids(2))
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if
        ierr = nf90_inq_varid(file_id, "c", var_ids(1))
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_var(file_id, var_ids(1), c, start = (/ 1,1,PT /))
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_var(file_id, var_ids(2), F_tot, start = (/PT/))
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if


        !print *, "Sucessfully written netcdf file"

    end subroutine write_netcdf_parts

    subroutine close_netcdf(file_id)
        integer, intent(in) :: file_id
        call check(nf90_close(file_id))
    end subroutine close_netcdf


end module io
