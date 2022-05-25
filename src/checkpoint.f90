module checkpointing

    use netcdf
    use iso_fortran_env
    use grid
    use comms

    implicit none

contains

    subroutine write_checkpoint_file(arr3dim, arr2dim,temp2dim, arr1dim,prob, coeffs, cpo, &
                                     initial_conc, nx, ny, m1, m2, k, bfe, Cint, t, time_step, current_iter, &
                                     df_tol, random_seed, ierr)


        real(kind=real64), intent(IN), dimension(:, :, :) :: arr3dim
        real(kind=real64), intent(IN), dimension(:, :) :: arr2dim,temp2dim
        real(kind=real64), intent(IN), dimension(:) :: arr1dim
        real(kind=real64), intent(IN), dimension(:) :: coeffs
        integer, intent(in) :: nx, ny, cint, random_seed, current_iter
        real(kind=real64), intent(in) :: initial_conc, m1, m2, k, bfe, t, &
                                         time_step, df_tol
        !TYPE(run_data_type), INTENT(IN) :: run_data
        character(len=*), intent(in) :: cpo,prob
        character(LEN=1), dimension(9) :: dims = (/"x", "y", "t", "c", "a", "b", "f","n","m"/)
        character(LEN=12), dimension(15) :: comp_names = (/"initial_conc", &
                                                            "nx          ", "ny          ", "m1          ", &
                                                           "m2          ", "k           ", "bfe         ", "cint        " &
                                                           , "max_t       ", "time_step   ", "current_time", "df_tol      " &
                                                           , "random_seed ", "cpo         ","problem     "/)
        integer :: file_id, i
        integer, intent(out) :: ierr
        integer :: ndims = 9
        integer, dimension(9) :: sizes, dim_ids
        integer, dimension(5) :: var_ids

        ! Acquiring the size of the dimensions
        sizes(1:3) = shape(arr3dim)
        sizes(4) = size(coeffs)
        sizes(5:6) = shape(arr2dim)
        sizes(7) = size(arr1dim)
        sizes(8:9) = shape(temp2dim)

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

        ierr = nf90_def_var(file_id, "mu", NF90_REAL, dim_ids(5:6), var_ids(3))
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_def_var(file_id, "ftot", NF90_REAL, dim_ids(7), var_ids(4))
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_def_var(file_id, "tempg", NF90_REAL, dim_ids(8:9), var_ids(5))
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

        ierr = nf90_put_var(file_id, var_ids(3), arr2dim)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_var(file_id, var_ids(4), arr1dim)
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if

        ierr = nf90_put_var(file_id, var_ids(5), temp2dim)
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

    subroutine read_checkpoint_metadata(fn, prob,initial_conc, &
                                   coeffs, Nx, Ny, M1, M2, k, bfe, Cint, cpo, t, delta_t, df_tol, &
                                  current_iter, random_seed, use_input, ierr)

        integer, intent(inout) :: nx, ny, cint, random_seed
        integer, intent(in) :: use_input
        character(len=*), intent(inout) :: cpo,prob
        character(len=*), intent(in) :: fn
        real(kind=real64), intent(inout) :: initial_conc, m1, m2, k, bfe, &
                                            t, delta_t, df_tol
        integer, intent(inout) :: current_iter
        real(kind=real64), intent(inout), dimension(:), allocatable :: coeffs
        character(LEN=1), dimension(9) :: dims = (/"x", "y", "t", "c", "a", "b", "f","n","m"/)
        character(LEN=12), dimension(15) :: comp_names = (/"initial_conc", &
                                                           "nx          ", "ny          ", "m1          ", &
                                                           "m2          ", "k           ", "bfe         ", "cint        " &
                                                           , "max_t       ", "time_step   ", "current_time", "df_tol      " &
                                                           , "random_seed ", "cpo         ","problem     "/)
        integer :: file_id, i
        integer, intent(out) :: ierr
        integer :: ndims = 9
        integer, dimension(9) :: sizes, dim_ids
        integer, dimension(5) :: var_ids

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

    subroutine read_checkpoint_data(arr3dim, arr2dim,temp2dim, arr1dim, fn,use_input, ierr)

        integer, intent(in) :: use_input
        character(len=*), intent(in) :: fn
        real(kind=real64), intent(inout), dimension(:, :, :), allocatable :: arr3dim
        real(kind=real64), intent(inout), dimension(:, :), allocatable :: arr2dim
        real(kind=real64), intent(inout), dimension(:, :), allocatable :: temp2dim
        real(kind=real64), intent(inout), dimension(:), allocatable :: arr1dim
        character(LEN=1), dimension(9) :: dims = (/"x", "y", "t", "c", "a", "b", "f","n","m"/)
        character(LEN=12), dimension(15) :: comp_names = (/"initial_conc", &
                                                           "nx          ", "ny          ", "m1          ", &
                                                           "m2          ", "k           ", "bfe         ", "cint        " &
                                                           , "max_t       ", "time_step   ", "current_time", "df_tol      " &
                                                           , "random_seed ", "cpo         ","problem     "/)
        integer :: file_id, i
        integer, intent(out) :: ierr
        integer :: ndims = 9
        integer, dimension(9) :: sizes, dim_ids
        integer, dimension(5) :: var_ids

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
                start=(/grid_domain_start(1),grid_domain_start(2),2/))
        if (ierr /= nf90_noerr) then
            print *, trim(nf90_strerror(ierr))
            return
        end if


        ! get var ids and read data
        if(my_rank == 0) then

            ierr = nf90_inq_varid(file_id, "mu", var_ids(3))
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_var(file_id, var_ids(3), arr2dim)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_inq_varid(file_id, "ftot", var_ids(4))
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            ierr = nf90_get_var(file_id, var_ids(4), arr1dim)
            if (ierr /= nf90_noerr) then
                print *, trim(nf90_strerror(ierr))
                return
            end if

            if (use_input == 0) then
                ierr = nf90_inq_varid(file_id, "tempg", var_ids(5))
                if (ierr /= nf90_noerr) then
                    print *, trim(nf90_strerror(ierr))
                    return
                end if

                ierr = nf90_get_var(file_id, var_ids(5), temp2dim)
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
