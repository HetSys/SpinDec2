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
      integer, intent (in) :: status
      
      if(status /= nf90_noerr) then 
        print *, trim(nf90_strerror(status))
        stop ""
      end if
    end subroutine check  
    
    !************************************************************************
    !> write_netcdf
    !!
    !! subroutine to perform writing ouput file (netcdf)
    !!
    !!@param filename : filename to store data in
    !!@param c : order parameter field over all time
    !!@param F_tot : total free energy over time
    !!@param coeffs : polynomial coefficients for f(c)
    !!@param N_c : number of coefficients
    !!@params N_x, N_y, N_t : number of  discretizations  in x, y and time
    !!@param dt : forward Euler time-step
    !!@param c_0 : average concentration
    !!@params M_A, M_B : Atomic mobilities of species A and B
    !!@param kappa : gradient term coefficient
    !************************************************************************
    !
    subroutine write_netcdf(c, F_tot, coeffs, N_c, N_x, N_y, N_t, dt, c_0, M_A, M_B, kappa)

        integer, intent(in) :: N_x, N_y, N_t, N_c
        real(kind=real64), dimension(N_x, N_y, N_t), intent(in) :: c
        real(kind=real64), dimension(N_t), intent(in) :: F_tot
        real(kind=real64), dimension(N_c), intent(in) :: coeffs
        integer, dimension(3) :: size_c
        integer, dimension(3) :: c_dim_ids
        integer :: size_F_tot, size_coeffs
        integer :: F_tot_dim_ids, coeffs_dim_ids
        integer :: c_var_id, F_tot_var_id, coeffs_var_id
        character(len=*), dimension(3), parameter :: c_dims=(/"c_x", "c_y", "c_t"/)
        character(len=*), parameter :: F_tot_dims = "F_t", coeffs_dims = "C_i"
        character(len=*), parameter :: filename = 'Visualise/test_netcdf.nc'
        real(kind=real64), intent(in) :: dt, M_A, M_B, kappa, c_0
        integer :: k, file_id
        
    
        !define size of arrays
        size_c = shape(c)
        size_F_tot = size(F_tot)
        size_coeffs = size(coeffs)
      

        !Create the file, overwriting if it exists
        call check(nf90_create(filename, NF90_CLOBBER, file_id))
    
   
        !write in the dimensions for the variables
        do k = 1, 3

          call check(nf90_def_dim(file_id, c_dims(k), size_c(k), c_dim_ids(k)))
          
        end do

        call check(nf90_def_dim(file_id, F_tot_dims, size_F_tot, F_tot_dim_ids))

        call check(nf90_def_dim(file_id, coeffs_dims, size_coeffs, coeffs_dim_ids))
  
        print *, "Dimensions Written"


        !Define variable type, matching the arrays
        call check (nf90_def_var(file_id, "c", NF90_DOUBLE, c_dim_ids, c_var_id))
        
        call check (nf90_def_var(file_id, "F_tot", NF90_DOUBLE, F_tot_dim_ids, F_tot_var_id))

        call check (nf90_def_var(file_id, "coeffs", NF90_DOUBLE, coeffs_dim_ids, coeffs_var_id))
       
        print *, "Variable types defined"
  
        ! Setting up meta variables of simulation

        !Time step
        call check(nf90_put_att(file_id, NF90_GLOBAL, "dt", dt))
        !Number of time iterations
        call check(nf90_put_att(file_id, NF90_GLOBAL, "N_t", N_t))
        !Grid size (x)
        call check(nf90_put_att(file_id, NF90_GLOBAL, "N_x", N_x))
        !Grid size (y)
        call check(nf90_put_att(file_id, NF90_GLOBAL, "N_y", N_y))
        !Mobility
        call check(nf90_put_att(file_id, NF90_GLOBAL, "M_A", M_A))
        call check(nf90_put_att(file_id, NF90_GLOBAL, "M_B", M_B))
        !gradient term coefficient
        call check(nf90_put_att(file_id, NF90_GLOBAL, "kappa", kappa))
        !initial average concentration
        call check(nf90_put_att(file_id, NF90_GLOBAL, "c_0", c_0))
       
        !Finish defining metadata
        call check(nf90_enddef(file_id))
       
        ! Actually write the variables
        call check(nf90_put_var(file_id, c_var_id, c))
        call check(nf90_put_var(file_id, F_tot_var_id, F_tot))
        call check(nf90_put_var(file_id, coeffs_var_id, coeffs))
        
        !close the file
        call check(nf90_close(file_id))
  
        print*, "SUCESSFULLY WRITTEN NETCDF FILE"
       
    
    end subroutine write_netcdf
    
end module io
