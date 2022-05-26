module checkpointing

    use netcdf
    use iso_fortran_env
    use grid
    use comms

    implicit none

contains

    !> Subroutine to write a checkpoint file
    !! @param arr3dim The 3 dimentional array to write, usually concentration grid
    !! @param temp2dim The 2 dimensional array to write, usually the temperature grid
    !! @param arr1dim The 1 dimensional array to write, usually the free energy over time
    !! @param prob Problem that will be computed (spectral, constant, nontemp, temp)
    !! @param coeffs Coefficients of the polynomial used for bulk potential
    !! @param cpo Name of the checkpoint out file
    !! @param initial_conc Initial average concentration of the grid
    !! @param nx X dimension of the grid
    !! @param ny Y dimension of the grid
    !! @param m1 Mobility of the first species
    !! @param m2 Mobility of the second species
    !! @param k The charecteristic lengthscale of transition region
    !! @param bfe Bulk free energy paramter
    !! @param Cint Checkpoint interval
    !! @param t Final time to run simulation till
    !! @param time_step Time step
    !! @param current_iter The interation at which the checkpoint was written
    !! @param random_seed Seed used for the simulation
    !! @param lastw The interation where the output was last written
    !! @param write_freq The frequency at which frames are written
    !! @param write_freq_c A counter that is saved so the calculation can continue
    !! @param write_count A counter that is saved so the calculation can continue
    !! @param singl A flag that tells the code to save in single single precision
    !! @param ierr A flag that checks if an error occured in writing
    subroutine write_checkpoint_file(arr3dim,temp2dim, arr1dim,prob, coeffs, cpo, &
                                     initial_conc, nx, ny, m1, m2, k, bfe, Cint, t, time_step, current_iter, &
                                     df_tol, random_seed,lastw,write_freq,write_freq_c,singl,write_count, stab,ierr)


        real(kind=real64), intent(IN), dimension(:, :, :) :: arr3dim
        real(kind=real64), intent(IN), dimension(:, :) :: temp2dim
        real(kind=real64), intent(IN), dimension(:) :: arr1dim
        real(kind=real64), intent(IN), dimension(:) :: coeffs
        integer, intent(in) :: nx, ny, cint, random_seed, current_iter,lastw,&
                                write_freq,write_freq_c,singl,write_count
        real(kind=real64), intent(in) :: initial_conc, m1, m2, k, bfe, t, &
                                         time_step, df_tol,stab
        !TYPE(run_data_type), INTENT(IN) :: run_data
        character(len=*), intent(in) :: cpo,prob
        character(LEN=1), dimension(7) :: dims = (/"x", "y", "t", "c", "f","n","m"/)
        character(LEN=12), dimension(21) :: comp_names = (/"initial_conc", &
                                                           "nx          ", "ny          ", "m1          ", &
                                                           "m2          ", "k           ", "bfe         ", "cint        " &
                                                           , "max_t       ", "time_step   ", "current_time", "df_tol      " &
                                                           , "random_seed ", "cpo         ","problem     ","last_write  "&
                                                           ,"write_freq  ","write_freq_c","singl       ","write_count "&
                                                           ,"stab        "/)
        integer :: file_id, i
        integer, intent(out) :: ierr
        integer :: ndims = 7
        integer, dimension(7) :: sizes, dim_ids
        integer, dimension(4) :: var_ids

        ! Acquiring the size of the dimensions
        sizes(1:3) = shape(arr3dim)
        sizes(4) = size(coeffs)
        sizes(5) = size(arr1dim)
        sizes(6:7) = shape(temp2dim)

        ! Opening a file
        ierr = nf90_create(cpo ,IOR(NF90_NETCDF4,NF90_CLOBBER), file_id)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ! Defining dimensions
        do i = 1, ndims
            ierr = nf90_def_dim(file_id, dims(i), sizes(i), dim_ids(i))
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if
        end do

        !Defining variables
        ierr = nf90_def_var(file_id, "data", NF90_REAL, dim_ids(1:3), var_ids(1))
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_def_var(file_id, "coeffs", NF90_REAL, dim_ids(4), var_ids(2))
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_def_var(file_id, "ftot", NF90_REAL, dim_ids(5), var_ids(3))
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_def_var(file_id, "tempg", NF90_REAL, dim_ids(6:7), var_ids(4))
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ! Adding global attributes (run data)
        ierr = nf90_put_att(file_id, NF90_GLOBAL, trim(comp_names(1)), initial_conc)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_att(file_id, NF90_GLOBAL, trim(comp_names(2)), nx)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_att(file_id, NF90_GLOBAL, trim(comp_names(3)), ny)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_att(file_id, NF90_GLOBAL, trim(comp_names(4)), m1)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_att(file_id, NF90_GLOBAL, trim(comp_names(5)), m2)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_att(file_id, NF90_GLOBAL, trim(comp_names(6)), k)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_att(file_id, NF90_GLOBAL, trim(comp_names(7)), bfe)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_att(file_id, NF90_GLOBAL, trim(comp_names(8)), cint)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_att(file_id, NF90_GLOBAL, trim(comp_names(9)), t)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_att(file_id, NF90_GLOBAL, trim(comp_names(10)), time_step)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_att(file_id, NF90_GLOBAL, trim(comp_names(11)), current_iter)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_att(file_id, NF90_GLOBAL, trim(comp_names(12)), df_tol)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_att(file_id, NF90_GLOBAL, trim(comp_names(13)), random_seed)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_att(file_id, NF90_GLOBAL, trim(comp_names(14)), cpo)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_att(file_id, NF90_GLOBAL, trim(comp_names(15)), prob)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_att(file_id, NF90_GLOBAL, trim(comp_names(16)), lastw)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_att(file_id, NF90_GLOBAL, trim(comp_names(17)), write_freq)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_att(file_id, NF90_GLOBAL, trim(comp_names(18)), write_freq_c)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_att(file_id, NF90_GLOBAL, trim(comp_names(19)), singl)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_att(file_id, NF90_GLOBAL, trim(comp_names(20)), write_count)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_att(file_id, NF90_GLOBAL, trim(comp_names(21)), stab)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ! Ending definition mode
        ierr = nf90_enddef(file_id)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ! Writing the data
        ierr = nf90_put_var(file_id, var_ids(1), arr3dim)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_var(file_id, var_ids(2), coeffs)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_var(file_id, var_ids(3), arr1dim)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_var(file_id, var_ids(4), temp2dim)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ! Closeing the file
        ierr = nf90_close(file_id)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        print *, "Writing checkpoint successful at iter", current_iter

    end subroutine write_checkpoint_file


    !> Subroutine to read metadata from a checkpoint file
    !! @param prob Problem that will be computed (spectral, constant, nontemp, temp)
    !! @param coeffs Coefficients of the polynomial used for bulk potential
    !! @param cpo Name of the checkpoint out file
    !! @param initial_conc Initial average concentration of the grid
    !! @param nx X dimension of the grid
    !! @param ny Y dimension of the grid
    !! @param m1 Mobility of the first species
    !! @param m2 Mobility of the second species
    !! @param k The charecteristic lengthscale of transition region
    !! @param bfe Bulk free energy paramter
    !! @param Cint Checkpoint interval
    !! @param t Final time to run simulation till
    !! @param time_step Time step
    !! @param current_iter The interation at which the checkpoint was written
    !! @param random_seed Seed used for the simulation
    !! @param lastw The interation where the output was last written
    !! @param write_freq The frequency at which frames are written
    !! @param write_freq_c A counter that is saved so the calculation can continue
    !! @param write_count A counter that is saved so the calculation can continue
    !! @param singl A flag that tells the code to save in single single precision
    !! @param fn Name of the file to read from
    !! @param use_input This is a flag that indicated is the input should override the checkpoint metadata
    !! @param ierr A flag that checks if an error occured in reading
    !! @param stab A stabilization term that allows for bigger timesteps at the const of inital accuracy (spectral)
    subroutine read_checkpoint_metadata(fn, prob,initial_conc, &
                                   coeffs, Nx, Ny, M1, M2, k, bfe, Cint, cpo, t, delta_t, df_tol, &
                                  current_iter, random_seed, use_input,&
                                  lastw,write_freq,write_freq_c,singl,write_count,stab, ierr)

        integer, intent(inout) :: nx, ny, cint, random_seed,write_freq,write_freq_c,singl,write_count
        integer, intent(in) :: use_input
        character(len=*), intent(inout) :: cpo,prob
        character(len=*), intent(in) :: fn
        real(kind=real64), intent(inout) :: initial_conc, m1, m2, k, bfe, &
                                            t, delta_t, df_tol,stab
        integer, intent(inout) :: current_iter,lastw
        real(kind=real64), intent(inout), dimension(:), allocatable :: coeffs
        character(LEN=1), dimension(7) :: dims = (/"x", "y", "t", "c", "f","n","m"/)
        character(LEN=12), dimension(21) :: comp_names = (/"initial_conc", &
                                                           "nx          ", "ny          ", "m1          ", &
                                                           "m2          ", "k           ", "bfe         ", "cint        " &
                                                           , "max_t       ", "time_step   ", "current_time", "df_tol      " &
                                                           , "random_seed ", "cpo         ","problem     ","last_write  "&
                                                           ,"write_freq  ","write_freq_c","singl       ","write_count "&
                                                           ,"stab        "/)
        integer :: file_id, i
        integer, intent(out) :: ierr
        integer :: ndims = 7
        integer, dimension(7) :: sizes, dim_ids
        integer, dimension(4) :: var_ids

        ! Open file
        ierr = nf90_open(fn, NF90_NOWRITE, file_id)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ! Get dimensions
        do i = 1, ndims
            ierr = nf90_inq_dimid(file_id, dims(i), dim_ids(i))
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if
        end do

        ! Get dimension lengths
        do i = 1, ndims
            ierr = nf90_inquire_dimension(file_id, dim_ids(i), len=sizes(i))
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if
        end do

        if (use_input == 0) then
            deallocate (coeffs)
            allocate (coeffs(sizes(4)))
            ierr = nf90_inq_varid(file_id, "coeffs", var_ids(2))
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_var(file_id, var_ids(2), coeffs)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_att(file_id, NF90_GLOBAL, trim(comp_names(1)), initial_conc)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_att(file_id, NF90_GLOBAL, trim(comp_names(2)), nx)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_att(file_id, NF90_GLOBAL, trim(comp_names(3)), ny)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_att(file_id, NF90_GLOBAL, trim(comp_names(4)), m1)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_att(file_id, NF90_GLOBAL, trim(comp_names(5)), m2)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_att(file_id, NF90_GLOBAL, trim(comp_names(6)), k)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_att(file_id, NF90_GLOBAL, trim(comp_names(7)), bfe)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_att(file_id, NF90_GLOBAL, trim(comp_names(8)), cint)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_att(file_id, NF90_GLOBAL, trim(comp_names(9)), t)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_att(file_id, NF90_GLOBAL, trim(comp_names(10)), delta_t)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_att(file_id, NF90_GLOBAL, trim(comp_names(12)), df_tol)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_att(file_id, NF90_GLOBAL, trim(comp_names(13)), random_seed)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_att(file_id, NF90_GLOBAL, trim(comp_names(14)), cpo)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_att(file_id, NF90_GLOBAL, trim(comp_names(15)), prob)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_att(file_id, NF90_GLOBAL, trim(comp_names(16)), lastw)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_att(file_id, NF90_GLOBAL, trim(comp_names(17)), write_freq)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_att(file_id, NF90_GLOBAL, trim(comp_names(18)), write_freq_c)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_att(file_id, NF90_GLOBAL, trim(comp_names(19)), singl)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_att(file_id, NF90_GLOBAL, trim(comp_names(20)), write_count)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_att(file_id, NF90_GLOBAL, trim(comp_names(21)), stab)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

        end if

        ierr = nf90_get_att(file_id, NF90_GLOBAL, trim(comp_names(11)), current_iter)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        !print*, current_time

        ierr = nf90_close(file_id)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        print *, "Reading checkpoint metadata successful"

    end subroutine read_checkpoint_metadata
    !> Subroutine to read data from a checkpoint file
    !! @param arr3dim The 3 dimentional array to write, usually concentration grid
    !! @param temp2dim The 2 dimensional array to write, usually the temperature grid
    !! @param arr1dim The 1 dimensional array to write, usually the free energy over time
    !! @param fn Name of the file to read from
    !! @param use_input This is a flag that indicated is the input should override the checkpoint metadata
    !! @param ierr A flag that checks if an error occured in reading
    subroutine read_checkpoint_data(arr3dim,temp2dim, arr1dim, fn,use_input, ierr)

        integer, intent(in) :: use_input
        character(len=*), intent(in) :: fn
        real(kind=real64), intent(inout), dimension(:, :, :), allocatable :: arr3dim
        real(kind=real64), intent(inout), dimension(:, :), allocatable :: temp2dim
        real(kind=real64), intent(inout), dimension(:), allocatable :: arr1dim
        character(LEN=1), dimension(7) :: dims = (/"x", "y", "t", "c", "f","n","m"/)

        integer :: file_id, i
        integer, intent(out) :: ierr
        integer :: ndims = 7
        integer, dimension(7) :: sizes, dim_ids
        integer, dimension(4) :: var_ids

        ! Open file
        ierr = nf90_open(fn, NF90_NOWRITE, file_id)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ! Get dimensions
        do i = 1, ndims
            ierr = nf90_inq_dimid(file_id, dims(i), dim_ids(i))
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if
        end do

        ! Get dimension lengths
        do i = 1, ndims
            ierr = nf90_inquire_dimension(file_id, dim_ids(i), len=sizes(i))
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if
        end do

        ierr = nf90_inq_varid(file_id, "data", var_ids(1))
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if
        ierr = nf90_get_var(file_id, var_ids(1), local_grid_conc,&
                start=(/grid_domain_start(2),grid_domain_start(1),1/))
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if


        ! get var ids and read data
        if(my_rank == 0) then
            ierr = nf90_inq_varid(file_id, "ftot", var_ids(3))
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_var(file_id, var_ids(3), arr1dim)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            if (use_input == 0) then
                ierr = nf90_inq_varid(file_id, "tempg", var_ids(4))
                if (ierr /= nf90_noerr) then
                    print *, trim(nf90_strerror(ierr))
                    return
                end if

                ierr = nf90_get_var(file_id, var_ids(4), temp2dim)
                if (ierr /= nf90_noerr) then
                    print *, trim(nf90_strerror(ierr))
                    return
                end if
            end if
        end if

        !print*, current_time

        ierr = nf90_close(file_id)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        print *, "Reading checkpoint data successful"

    end subroutine read_checkpoint_data

end module checkpointing
