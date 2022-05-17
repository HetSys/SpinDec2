! Code adapted from PX425 - Assignment 4 (2021) By Dr. David Quigley
module get_global_grid

    use iso_fortran_env
    
    implicit none

    contains
    
    ! Routine to collect all contributions to the global grid onto      !
    ! rank zero for visualisation.                                      !
    subroutine comms_get_global_grid()
       
        ! comms buffer
        real(real64),allocatable,dimension(:) :: combuff
    
        ! MPI Status
        integer,dimension(MPI_STATUS_SIZE) :: status
    
        ! Information on the remote domain
        integer,dimension(2) :: remote_domain_start= (/1,1/)
    
        ! Loop counters and error flag
        integer :: ix,iy,ixg,iyg,ierr,ip
    
        ! Just use local grid if running on one processor
        if (p==1) then
            global_grid_conc=local_grid_conc
            return
        end if
    
        if (my_rank==0) then
    
            ! Rank 0 first fills out his part of the global grid
            do iy = 1,grid_domain_size
                do ix = 1,grid_domain_size
    
                   ! Global indices
                    ixg = ix+grid_domain_start(x)-1
                    iyg = iy+grid_domain_start(y)-1
    
                    global_grid_conc(ixg,iyg) = local_grid_conc(ix,iy)
    
                end do
            end do
    
        end if
        
        ! Allocate buffer
        allocate(combuff(1:grid_domain_size),stat=ierr)
        if (ierr/=0) stop 'Error allocating combuff in comms_get_global_grid'


        if (my_rank==0) then

            ! Now loop over all other ranks receiving their data
            do ip = 1,p-1

            ! First receive remote_domain_start from rank ip
            ! Insert appropriate MPI call here
                call MPI_Recv(remote_domain_start,2,MPI_INT,ip,999,MPI_COMM_WORLD,status,ierr)

            ! Loop over columns within a domain
                do iy = 1,grid_domain_size

                    ! Receive this column from rank ip
                    ! Insert appropriate MPI call here
                    call MPI_Recv(combuff,grid_domain_size,MPI_DOUBLE_PRECISION,ip,888,MPI_COMM_WORLD,status,ierr)

                    do ix = 1,grid_domain_size

                        ! Global indices
                        ixg = ix+remote_domain_start(x)-1
                        iyg = iy+remote_domain_start(y)-1

                        ! Store in global_grid_conc
                        global_grid_conc(ixg,iyg) = combuff(ix)

                    end do ! elements in column

                end do ! columns

            end do ! processors

        else ! not rank 0

       ! All other processors must send the data rank 0 needs

       ! Send grid_domain_start to rank 0
       ! Insert appropriate MPI call here
            call MPI_Send(grid_domain_start,2,MPI_INT,0,999,MPI_COMM_WORLD,ierr)

            do iy = 1,grid_domain_size

        ! Insert appropriate MPI call here
                call MPI_Send(local_grid_conc(:,iy),grid_domain_size,MPI_DOUBLE_PRECISION,0,888,MPI_COMM_WORLD,ierr)

            end do

        end if

       ! Free memory
       deallocate(combuff,stat=ierr)
       if (ierr/=0) stop "Error deallocating combuff in comms_get_global_grid!"
 
end module get_global_grid
