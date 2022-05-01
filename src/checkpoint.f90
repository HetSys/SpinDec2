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
        SUBROUTINE write_checkpoint_file(arr3dim, coeffs,cpo,initial_conc,conc_std&
            ,nx,ny,m1,m2,k,bfe,Cint,t,time_step,current_time,df_tol,random_seed,ierr)
            real(kind=real64), INTENT(IN), DIMENSION(:,:,:) :: arr3dim
            real(kind=real64), INTENT(IN), DIMENSION(:) :: coeffs
            integer , intent(in) :: nx,ny,cint,random_seed
            real(kind=REAL64), intent(in) :: initial_conc, conc_std,m1,m2,k,bfe,t,time_step&
                ,df_tol,current_time
            !TYPE(run_data_type), INTENT(IN) :: run_data
            character(len = *),intent(in) :: cpo
            CHARACTER(LEN=1), DIMENSION(4) :: dims = (/"x", "y", "t","c"/)
            CHARACTER(LEN=12), DIMENSION(15) :: comp_names = (/"initial_conc",&
            "conc_std    ", "nx          ","ny          ","m1          ",&
            "m2          ","k           ","bfe         ","cint        "&
            ,"max_t       ","time_step   ","current_time","df_tol      "&
            ,"random_seed ","cpo         "/)
            INTEGER :: file_id, i
            INTEGER,intent(out) :: ierr
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
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(1)),initial_conc)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(2)),conc_std)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(3)),nx)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(4)),ny)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(5)),m1)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(6)),m2)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(7)),k)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(8)),bfe)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(9)),cint)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(10)),t)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(11)),time_step)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(12)),current_time)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(13)),df_tol)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF
            ierr = nf90_put_att(file_id, NF90_GLOBAL, TRIM(comp_names(14)),random_seed)
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


        SUBROUTINE read_checkpoint_in(arr3dim, fn, initial_conc,conc_std,coeffs&
            ,Nx,Ny,M1,M2,k,bfe,Cint,cpo,t,delta_t,df_tol,current_time,random_seed,use_input,ierr)
            integer , intent(inout) :: nx,ny,cint,random_seed
            integer,intent(in) :: use_input
            character(len=128),intent(inout) :: cpo
            character(len=*),intent(in) :: fn
            real(kind=REAL64), intent(inout) :: initial_conc, conc_std,m1,m2,k,bfe,t,delta_t&
                ,df_tol,current_time
            real(kind=real64), intent(inout),dimension(:),allocatable :: coeffs
            real(kind=real64), intent(out),dimension(:,:,:),allocatable :: arr3dim
            CHARACTER(LEN=1), DIMENSION(4) :: dims = (/"x", "y", "t","c"/)
            CHARACTER(LEN=12), DIMENSION(15) :: comp_names = (/"initial_conc",&
            "conc_std    ", "nx          ","ny          ","m1          ",&
            "m2          ","k           ","bfe         ","cint        "&
            ,"max_t       ","time_step   ","current_time","df_tol      "&
            ,"random_seed ","cpo         "/)
            INTEGER ::  file_id, i
            INTEGER,intent(out) :: ierr
            INTEGER :: ndims = 4
            INTEGER, DIMENSION(4) :: sizes, dim_ids
            INTEGER, DIMENSION(2) :: var_ids

            !Open file
            ierr = nf90_open(fn, NF90_NOWRITE, file_id)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF

            !Get dimensions
            do i = 1,ndims
                ierr = nf90_inq_dimid(file_id,dims(i),dim_ids(i))
                IF (ierr /= nf90_noerr) THEN
                    PRINT*, TRIM(nf90_strerror(ierr))
                    RETURN
                END IF
            end do

            !Get dimension lengths
            do i = 1,ndims
                ierr = nf90_inquire_dimension(file_id, dim_ids(i), len = sizes(i))
                IF (ierr /= nf90_noerr) THEN
                    PRINT*, TRIM(nf90_strerror(ierr))
                    RETURN
                END IF
            end do

            !allocate
            allocate(arr3dim(sizes(1),sizes(2),sizes(3)))


            !get var ids and read data
            ierr = nf90_inq_varid(file_id,"data",var_ids(1))
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF

            ierr = nf90_get_var(file_id,var_ids(1),arr3dim)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF


            if(use_input == 0) then
                deallocate(coeffs)
                allocate(coeffs(sizes(4)))
                ierr = nf90_inq_varid(file_id,"coeffs",var_ids(2))
                IF (ierr /= nf90_noerr) THEN
                    PRINT*, TRIM(nf90_strerror(ierr))
                    RETURN
                END IF

                ierr = nf90_get_var(file_id,var_ids(2),coeffs)
                IF (ierr /= nf90_noerr) THEN
                    PRINT*, TRIM(nf90_strerror(ierr))
                    RETURN
                END IF


                ierr = nf90_get_att(file_id, NF90_GLOBAL, TRIM(comp_names(1)),initial_conc)
                IF (ierr /= nf90_noerr) THEN
                    PRINT*, TRIM(nf90_strerror(ierr))
                    RETURN
                END IF
                ierr = nf90_get_att(file_id, NF90_GLOBAL, TRIM(comp_names(2)),conc_std)
                IF (ierr /= nf90_noerr) THEN
                    PRINT*, TRIM(nf90_strerror(ierr))
                    RETURN
                END IF
                ierr = nf90_get_att(file_id, NF90_GLOBAL, TRIM(comp_names(3)),nx)
                IF (ierr /= nf90_noerr) THEN
                    PRINT*, TRIM(nf90_strerror(ierr))
                    RETURN
                END IF
                ierr = nf90_get_att(file_id, NF90_GLOBAL, TRIM(comp_names(4)),ny)
                IF (ierr /= nf90_noerr) THEN
                    PRINT*, TRIM(nf90_strerror(ierr))
                    RETURN
                END IF
                ierr = nf90_get_att(file_id, NF90_GLOBAL, TRIM(comp_names(5)),m1)
                IF (ierr /= nf90_noerr) THEN
                    PRINT*, TRIM(nf90_strerror(ierr))
                    RETURN
                END IF
                ierr = nf90_get_att(file_id, NF90_GLOBAL, TRIM(comp_names(6)),m2)
                IF (ierr /= nf90_noerr) THEN
                    PRINT*, TRIM(nf90_strerror(ierr))
                    RETURN
                END IF
                ierr = nf90_get_att(file_id, NF90_GLOBAL, TRIM(comp_names(7)),k)
                IF (ierr /= nf90_noerr) THEN
                    PRINT*, TRIM(nf90_strerror(ierr))
                    RETURN
                END IF
                ierr = nf90_get_att(file_id, NF90_GLOBAL, TRIM(comp_names(8)),bfe)
                IF (ierr /= nf90_noerr) THEN
                    PRINT*, TRIM(nf90_strerror(ierr))
                    RETURN
                END IF
                ierr = nf90_get_att(file_id, NF90_GLOBAL, TRIM(comp_names(9)),cint)
                IF (ierr /= nf90_noerr) THEN
                    PRINT*, TRIM(nf90_strerror(ierr))
                    RETURN
                END IF
                ierr = nf90_get_att(file_id, NF90_GLOBAL, TRIM(comp_names(10)),t)
                IF (ierr /= nf90_noerr) THEN
                    PRINT*, TRIM(nf90_strerror(ierr))
                    RETURN
                END IF
                ierr = nf90_get_att(file_id, NF90_GLOBAL, TRIM(comp_names(11)),delta_t)
                IF (ierr /= nf90_noerr) THEN
                    PRINT*, TRIM(nf90_strerror(ierr))
                    RETURN
                END IF
                ierr = nf90_get_att(file_id, NF90_GLOBAL, TRIM(comp_names(13)),df_tol)
                IF (ierr /= nf90_noerr) THEN
                    PRINT*, TRIM(nf90_strerror(ierr))
                    RETURN
                END IF
                ierr = nf90_get_att(file_id, NF90_GLOBAL, TRIM(comp_names(14)),random_seed)
                IF (ierr /= nf90_noerr) THEN
                    PRINT*, TRIM(nf90_strerror(ierr))
                    RETURN
                END IF
                ierr = nf90_get_att(file_id, NF90_GLOBAL, TRIM(comp_names(15)),cpo)
                IF (ierr /= nf90_noerr) THEN
                    PRINT*, TRIM(nf90_strerror(ierr))
                    RETURN
                END IF
            end if

            ierr = nf90_get_att(file_id, NF90_GLOBAL, TRIM(comp_names(12)),current_time)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF

            !print*, current_time



            ierr = nf90_close(file_id)
            IF (ierr /= nf90_noerr) THEN
                PRINT*, TRIM(nf90_strerror(ierr))
                RETURN
            END IF

            PRINT*, "Reading checkpoint successful"


        END subroutine read_checkpoint_in



    END MODULE checkpointing
