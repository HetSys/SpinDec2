!> Comms module contains all routines which interact
!! with MPI libraries. Derived from PX425 Assignment 3
!! of Term 2 in 2021-2022
!! Original code created by D.Quigley - November 2012      
module comms

    use iso_fortran_env, dp => real64
    use mpi
    use grid

    implicit none

    private

    public:: comms_initialise
    public:: comms_processor_map
    public:: comms_finalise
    public:: comms_halo_swaps
    public:: comms_get_global_F
    public:: comms_get_global_grid

    public:: p
    public:: my_rank
    public:: my_rank_coords
    public:: my_rank_neighbours

    integer:: p
    integer:: my_rank
    integer, dimension(2):: my_rank_coords
    integer, dimension(4):: my_rank_neighbours
    integer:: cart_comm  ! cartesian communicator

contains

    !> Subroutine to initialise MPI, get the communicator size p 
    !! and the rank my_rank of the current process within 
    !! that communicator.
    subroutine comms_initialise()

        integer:: ierr,prov

        !call mpi_init(ierr)
        call mpi_init_thread(MPI_THREAD_FUNNELED,prov,ierr)

        call mpi_comm_rank(mpi_comm_world, my_rank, ierr)
        call mpi_comm_size(mpi_comm_world, p, ierr)

    end subroutine comms_initialise

    !> Subroutine to map our p processors into a 2D Cartesian grid of    
    !! dimension proot by proot where proot = sqrt(p).                                                                                      !
    !! Should populate the arrays my_rank_cooords, which contains the    
    !! location of the current MPI task within the processor grid, and   
    !! my_rank_neighbours, which contains (in the order left, right,     
    !! down and up) ranks of neighbouring MPI tasks on the grid with     
    !! which the current task will need to communicate.                  
    subroutine comms_processor_map()

        ! setting up cartesian communicator

        integer, parameter:: ndims = 2
        logical:: reorder = .false.
        logical, dimension(2):: pbc = (/.true., .true./)
        integer, dimension(2):: dims

        integer:: proot
        integer:: ierr

        proot = int(real(sqrt(real(p, kind = dp)), kind = dp) + 0.5)

        dims(:) = proot

        ! create the cartesian communicator
        call mpi_cart_create(mpi_comm_world, ndims, dims, pbc, reorder, cart_comm, ierr)

        ! get rank within the new communicator
        call mpi_comm_rank(cart_comm, my_rank, ierr)

        ! get and store the coordinates of the current mpi task
        call mpi_cart_coords(cart_comm, my_rank, ndims, my_rank_coords, ierr)

        ! rank of neighbouring tasks in the four directions
        call mpi_cart_shift(cart_comm, 0, 1, my_rank_neighbours(4), my_rank_neighbours(3), ierr)

        call mpi_cart_shift(cart_comm, 1, 1, my_rank_neighbours(1), my_rank_neighbours(2), ierr)

    end subroutine comms_processor_map

    !> Subroutine to compute the global free energy of the grid by     
    !! summing over all values of loca free energy, and storing the     
    !! result in global_F.
    !! @param local_F Local free energy of current rank (input)
    !! @param global F Global free energy (output)                                                   
    subroutine comms_get_global_F(local_F,global_F)

        real(kind=dp),intent(in)  :: local_F
        real(kind=dp),intent(out) :: global_F

        integer :: ierr              ! Error flag

        call MPI_Reduce(local_F,global_F,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    end subroutine comms_get_global_F

    !> Subroutine to send boundary spins on each side of the local grid  
    !! to neighbour processors, and to receive from those processors the 
    !! halo information needed to perform computations involving spins   
    !! on the boundary processors grid.      
    !! @param grid 2D Grid of the current process
    !! @param grid_halo 2D Array to store data from neighbours
    subroutine comms_halo_swaps(grid,grid_halo)

        ! send and receive concs on each side of the grid to neighbour processors

        real(dp), dimension(:,:), intent(in) :: grid
        real(dp), dimension(:,:), intent(out) :: grid_halo
        real(dp), allocatable, dimension(:):: sendbuf, recvbuf
        integer, dimension(mpi_status_size):: status
        integer :: ierr

        ! exit if not in parallel
        if (p == 1) then
            grid_halo(:, right) = grid(1,:)
            grid_halo(:, left) = grid(grid_domain_size, :)
            grid_halo(:, up) = grid(:,grid_domain_size)
            grid_halo(:, down) = grid(:,1)
            return
        end if

        ! allocate buffers
        allocate (sendbuf(1:grid_domain_size), stat = ierr)
        if (ierr /= 0) then
            stop 'error allocating sendbuf x in comms_halo_swaps'
        end if

        allocate (recvbuf(1:grid_domain_size), stat = ierr)
        if (ierr /= 0) then
            stop 'error allocating recvbuf x in comms_halo_swaps'
        end if

        ! send left hand boundary elements of grid to my_rank_neighbours(left)
        ! and receive from my_rank_neighbours(right) into right of grid_halo
        sendbuf(:) = grid(1, :)
        call mpi_sendrecv(sendbuf, grid_domain_size, mpi_double_precision, my_rank_neighbours(left), 11, &
                          recvbuf, grid_domain_size, mpi_double_precision, my_rank_neighbours(right), 11, &
                          cart_comm, status, ierr)
        grid_halo(:, right) = recvbuf(:)

        ! send right hand boundary elements of grid to my_rank_neighbours(right)
        ! and receive from my_rank_neighbours(left) into left of grid_halo
        sendbuf(:) = grid(grid_domain_size, :)
        call mpi_sendrecv(sendbuf, grid_domain_size, mpi_double_precision, my_rank_neighbours(right), 12, &
                          recvbuf, grid_domain_size, mpi_double_precision, my_rank_neighbours(left), 12, &
                          cart_comm, status, ierr)
        grid_halo(:, left) = recvbuf(:)

        ! send bottom boundary elements of grid to my_rank_neighbours(down)
        ! and receive from my_rank_neighbours(up) into down of grid_halo
        sendbuf(:) = grid(:, grid_domain_size)
        call mpi_sendrecv(sendbuf, grid_domain_size, mpi_double_precision, my_rank_neighbours(down), 13, &
                          recvbuf, grid_domain_size, mpi_double_precision, my_rank_neighbours(up), 13, &
                          cart_comm, status, ierr)
        grid_halo(:, up) = recvbuf(:)

        ! send top boundary elements of grid to my_rank_neighbours(up)
        ! and receive from my_rank_neighbours(down) into up of grid_halo
        sendbuf(:) = grid(:, 1)
        call mpi_sendrecv(sendbuf, grid_domain_size, mpi_double_precision, my_rank_neighbours(up), 14, &
                          recvbuf, grid_domain_size, mpi_double_precision, my_rank_neighbours(down), 14, &
                          cart_comm, status, ierr)
        grid_halo(:, down) = recvbuf(:)

        deallocate (sendbuf, recvbuf, stat = ierr)
        if (ierr /= 0) then
            stop 'error releasing memory in comms_halo_swaps'
        end if

    end subroutine comms_halo_swaps

    !> Subroutine to finalise MPI.  
    subroutine comms_finalise()

        integer:: ierr

        call mpi_finalize(ierr)

    end subroutine comms_finalise

    !> Routine to collect all contributions to the global grid     
    !! onto rank zero.                                    
    subroutine comms_get_global_grid()

        ! comms buffer
        real(dp), allocatable, dimension(:):: combuff

        ! MPI Status
        integer, dimension(MPI_STATUS_SIZE):: status

        ! Information on the remote domain
        integer, dimension(2):: remote_domain_start = (/1, 1/)

        ! Loop counters and error flag
        integer:: ix, iy, ixg, iyg, ierr, ip

        ! Just use local grid if running on one processor
        if (p == 1) then
            global_grid_conc = local_grid_conc
            return
        end if


        if (my_rank == 0) then

            ! Rank 0 first fills out his part of the global grid
            do iy = 1, grid_domain_size
                do ix = 1, grid_domain_size

                    ! Global indices
                    ixg = ix+grid_domain_start(2) - 1
                    iyg = iy+grid_domain_start(1) - 1

                    global_grid_conc(ixg, iyg) = local_grid_conc(ix, iy)

                end do
            end do

        end if

        ! Allocate buffer
        allocate (combuff(1:grid_domain_size), stat = ierr)
        if (ierr /= 0) stop 'Error allocating combuff in comms_get_global_grid'

        if (my_rank == 0) then

            ! Now loop over all other ranks receiving their data
            do ip = 1, p-1

                ! First receive remote_domain_start from rank ip
                ! Insert appropriate MPI call here
                call MPI_Recv(remote_domain_start, 2, MPI_INT, ip, 999, MPI_COMM_WORLD, status, ierr)

                ! Loop over columns within a domain
                do iy = 1, grid_domain_size

                    ! Receive this column from rank ip
                    ! Insert appropriate MPI call here
                    call MPI_Recv(combuff, grid_domain_size, MPI_DOUBLE_PRECISION, ip, 888, &
                        MPI_COMM_WORLD, status, ierr)

                    do ix = 1, grid_domain_size

                        ! Global indices
                        ixg = ix+remote_domain_start(2) - 1
                        iyg = iy+remote_domain_start(1) - 1

                        ! Store in global_grid_conc
                        global_grid_conc(ixg, iyg) = combuff(ix)

                    end do  ! elements in column

                end do  ! columns

            end do  ! processors

        else  ! not rank 0

            ! All other processors must send the data rank 0 needs

            ! Send grid_domain_start to rank 0
            ! Insert appropriate MPI call here
            call MPI_Send(grid_domain_start, 2, MPI_INT, 0, 999, MPI_COMM_WORLD, ierr)

            do iy = 1, grid_domain_size

                ! Insert appropriate MPI call here
                call MPI_Send(local_grid_conc(:, iy), grid_domain_size, MPI_DOUBLE_PRECISION, 0, 888, &
                    MPI_COMM_WORLD, ierr)

            end do

        end if

        ! Free memory
        deallocate (combuff, stat = ierr)
        if (ierr /= 0) stop "Error deallocating combuff in comms_get_global_grid!"

    end subroutine comms_get_global_grid

end module commsc