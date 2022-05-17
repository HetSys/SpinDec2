! Code adapted from PX425 - Assignment 4 (2021) By Dr. David Quigley

module grid_mpi

    use iso_fortran_env
    use grid

    implicit none

    real(real64), dimension(:,:), allocatable :: global_grid_conc
    real(real64), dimension(:,:), allocatable :: local_grid_conc
    real(real64), dimension(:,:), allocatable :: halo_grid_conc
    
    integer :: grid_domain_size
    integer, dimension(2) :: grid_domain_start, grid_domain_end

    contains
    
    subroutine grid_initialise_local(Nx,Ny,p,my_rank,my_rank_coords)

        integer, intent(in) :: Nx, Ny, p
        integer, intent(in) :: my_rank
        integer, dimension(2), intent(in) :: my_rank_coords
        integer :: proot
        integer :: ierr

        proot = int(real(sqrt(real(p,kind=real64)),kind=real64)+0.5)

        grid_domain_size = (Nx*Ny)/proot
        grid_domain_start(:) = 1 + (my_rank_coords(:))*grid_domain_size
        grid_domain_end(:) = grid_domain_start(:) + grid_domain_size - 1

        if (my_rank == 0) then
            print*, "Size of each local processor grid is: ", grid_domain_size, " x ", grid_domain_size
        end if

        ! Set up four subdomains on current rank

        ! Allocate concnentration grid on current rank
        allocate(local_grid_conc(grid_domain_size,grid_domain_size),stat=ierr)
        if (ierr /= 0) then
            print*, "Error: allocating local_grid_con failed on rank ", my_rank
            stop
        end if

        ! Allocate four halo arrays
        allocate(halo_grid_conc(grid_domain_size,4),stat=ierr)
        if (ierr /= 0) then
            print*, "Error: allocating halo_grid_conc failed on rank ", my_rank
            stop
        end if

    end subroutine grid_initialise_local

    
    subroutine grid_initialise_global(Nx,Ny,c_min,c_max,my_rank)

        integer, intent(in) :: Nx, Ny, my_rank
        real(real64), intent(in) :: c_min, c_max
        integer :: ierr

        ! allocate grid array on rank zero
        if (my_rank == 0) then
            allocate(global_grid_conc(Nx,Ny),stat=ierr)
            if(ierr/=0) stop "Error: allocating global_grid_conc failed"

            ! Initiate grid using a uniform distribution
            call grid_init(global_grid_conc,Nx,Ny,c_min,c_max)
        end if

    end subroutine grid_initialise_global



    subroutine global_grid_deallocate(my_rank)

        integer, intent(in) :: my_rank
        integer :: ierr

        ! Deallocate memory allocated for global grid
        ! on rank 0
        if (my_rank == 0) then
            deallocate(global_grid_conc, stat=ierr)
            
            if (ierr /= 0) then
                print*, "Error: deallocating global_grid_conc failed on rank 0"
                stop
            end if
        end if

    end subroutine global_grid_deallocate


    subroutine local_grid_deallocate(my_rank)
        
        integer, intent(in) :: my_rank
        integer :: ierr

        ! Deallocate memory allocated for each local concentration grid
        deallocate(local_grid_conc, stat=ierr)

        if (ierr /= 0) then
            print*, "Error: deallocating local_grid_conc failed on rank ", my_rank
            stop
        end if

    end subroutine local_grid_deallocate


end module grid_mpi