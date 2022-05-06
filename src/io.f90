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
        character(len=*), parameter :: filename = 'CH_output.nc'
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

end module io
