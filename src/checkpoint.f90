MODULE checkpointing

    USE netcdf
    USE ISO_FORTRAN_ENV

    IMPLICIT NONE

    ! Creating a data type to hold run data
    TYPE run_data_type
        real(kind=real64) :: initial_conc
        real(kind=real64) :: conc_std
        integer :: nx
        integer :: ny
        real(kind=real64) :: m1
        real(kind=real64) :: m2
        real(kind=real64) :: k
        real(kind=real64) :: bfe
        integer :: Cint
        real(kind=real64) :: t
        real(kind=real64) :: time_step
        real(kind=real64) :: current_time
        real(kind=real64) :: df_tol
        integer :: random_seed
    END TYPE

    CONTAINS
        ! This subroutine was apdated from the example provided in the assignment brief
        SUBROUTINE write_checkpoint_file(arr3dim, coeffs,cpo,run_data)
            real(kind=real64), INTENT(IN), DIMENSION(:,:,:) :: arr3dim
            real(kind=real64), INTENT(IN), DIMENSION(:) :: coeffs
            TYPE(run_data_type), INTENT(IN) :: run_data
            character(len = *) :: cpo
            CHARACTER(LEN=1), DIMENSION(4) :: dims = (/"x", "y", "t","c"/)
            CHARACTER(LEN=12), DIMENSION(15) :: comp_names = (/"initial_conc",&
            "conc_std    ", "nx          ","ny          ","m1          ",&
            "m2          ","k           ","bfe         ","cint        "&
            ,"max_t       ","time_step   ","current_time","df_tol      "&
            ,"random_seed ","cpo         "/)
            INTEGER :: ierr, file_id, i
            INTEGER :: ndims = 4
            INTEGER, DIMENSION(4) :: sizes, dim_ids
            INTEGER, DIMENSION(2) :: var_ids

            ! Acquiring the size of the dimensions
            sizes(1:3) = SHAPE(arr3dim)
            sizes(4) = size(coeffs)


            ! Opening a file
            ierr = nf90_create(cpo, NF90_CLOBBER, file_id)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF

            ! Defining dimensions
            DO i = 1, ndims
                ierr = nf90_def_dim(file_id, dims(i), sizes(i), dim_ids(i))
                IF (ierr /= nf90_noerr) THEN
                    PRINT*, TRIM(nf90_strerror(ierr))
                    RETURN
                END IF
            END DO

            !Defining variables
            ierr = nf90_def_var(file_id, "data", NF90_REAL, dim_ids(1:3), var_ids(1))
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF

            ierr = nf90_def_var(file_id, "coeffs", NF90_REAL, dim_ids(4), var_ids(2))
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF

            ! Adding global attributes (run data)
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(1)),run_data%initial_conc)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(2)),run_data%conc_std)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(3)),run_data%nx)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(4)),run_data%ny)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(5)),run_data%m1)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(6)),run_data%m2)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(7)),run_data%k)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(8)),run_data%bfe)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(9)),run_data%cint)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(10)),run_data%t)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(11)),run_data%time_step)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(12)),run_data%current_time)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(13)),run_data%df_tol)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(14)),run_data%random_seed)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(15)),cpo)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF

            !Ending definition mode
            ierr = nf90_enddef(file_id)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF

            ! Writing the data
            ierr = nf90_put_var(file_id, var_ids(1), arr3dim)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF

            ierr = nf90_put_var(file_id, var_ids(2), coeffs)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF


            ! Closeing the file
            ierr = nf90_close(file_id)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF

            PRINT*, "Writing checkpoint successful"


        END SUBROUTINE write_checkpoint_file


    END MODULE checkpointing
