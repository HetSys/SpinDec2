    subroutine comms_halo_swaps()
        ! Send and receive concs on each side of the grid to neighbour processors

        integer, allocatable, dimension(:) :: sendbuf, recvbuf
        integer, dimension(mpi_status_size) :: status
        integer :: ix, iy, ierr

        ! Exit if not in parallel
        if (nprocs == 1) then
            return
        end if

        ! Allocate buffers
        allocate(sendbuf_x(1:Nx), stat=ierr)
        if (ierr /= 0) then
            stop 'Error allocating sendbuf x in comms_halo_swaps'
            ! TODO Remove after debugging
        end if

        allocate(recvbuf_x(1:Nx), stat=ierr)
        if (ierr /= 0) then
            stop 'Error allocating recvbuf x in comms_halo_swaps'
            ! TODO Remove after debugging
        end if

        allocate(sendbuf_y(1:Ny), stat=ierr)
        if (ierr /= 0) then
            stop 'Error allocating sendbuf y in comms_halo_swaps'
            ! TODO Remove after debugging
        end if

        allocate(recvbuf_y(1:Ny), stat=ierr)
        if (ierr /= 0) then
            stop 'Error allocating recvbuf y in comms_halo_swaps'
            ! TODO Remove after debugging
        end if

        ! Send left hand boundary elements of grid to my_rank_neighbours(left)
        ! and receive from my_rank_neighbours(right) into right of grid_halo
        sendbuf_y(:) = grid(1, :)
        call mpi_sendrecv(sendbuf_y, Ny, mpi_integer, my_rank_neighbours(left), 11, &
            recvbuf_y, Ny, mpi_integer, my_rank_neighbours(right), 11, cart_comm, status, ierr)
        grid_halo(:, right) = recvbuf_y(:)

        ! Send right hand boundary elements of grid to my_rank_neighbours(right)
        ! and receive from my_rank_neighbours(left) into left of grid_halo
        sendbuf_y(:) = grid(Ny, :)
        call mpi_sendrecv(sendbuf_y, Ny, mpi_integer, my_rank_neighbours(right), 12, &
            recvbuf_y, Ny, mpi_integer, my_rank_neighbours(left), 12, cart_comm, status, ierr)
        grid_halo(:, left) = recvbuf_y(:)

        ! Send bottom boundary elements of grid to my_rank_neighbours(down)
        ! and receive from my_rank_neighbours(up) into down of grid_halo
        sendbuf_x(:) = grid(:, Nx)
        call mpi_sendrecv(sendbuf_x, Nx, mpi_integer, my_rank_neighbours(down), 13, &
            recvbuf_x, Nx, mpi_integer, my_rank_neighbours(up), 13, cart_comm, status, ierr)
        grid_halo(:, down) = recvbuf_x(:)

        ! Send top boundary elements of grid to my_rank_neighbours(up)
        ! and receive from my_rank_neighbours(down) into up of grid_halo
        sendbuf_x(:) = grid_spin(:, 1)
        call mpi_sendrecv(sendbuf_x, Nx, mpi_integer, my_rank_neighbours(up), 14, &
            recvbuf_x, Nx, mpi_integer, my_rank_neighbours(down), 14, cart_comm, status, ierr)
        grid_halo(:, up) = recvbuf_x(:)

        deallocate(sendbuf_x, recvbuf_x, sendbuf_y, recvbuf_y, stat=ierr)
        if (ierr /= 0) then
            stop 'Error releasing memory in comms_halo_swaps'
            ! TODO Remove after debugging
        end if

    end subroutine comms_halo_swaps
