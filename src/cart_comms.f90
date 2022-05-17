! Code adapted from PX425 - Assignment 5 (2021) By Dr. David Quigley

module comms

    use iso_fortran_env
    use mpi
    ! use grid_mpi

    implicit none

    private

    public :: comms_initialise
    public :: comms_processor_map
    public :: comms_finalise

    public :: p
    public :: my_rank
    public :: my_rank_coords
    public :: my_rank_neighbours

    integer :: p
    integer :: my_rank
    integer, dimension(2) :: my_rank_coords
    integer, dimension(4) :: my_rank_neighbours
    integer :: cart_comm ! Cartesian communicator


    contains

    subroutine comms_initialise()

        integer :: proot
        integer :: ierr

        call MPI_Init(ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, p, ierr)

    end subroutine comms_initialise


    subroutine comms_processor_map()

        ! Setting up Cartesian communicator

        integer, parameter :: ndims = 2
        logical :: reorder = .FALSE.
        logical, dimension(2) :: pbc = (/.TRUE.,.TRUE./)
        integer, dimension(2) :: dims
        
        integer :: proot
        integer :: ierr

        proot = int(real(sqrt(real(p,kind=real64)),kind=real64)+0.5)

        dims(:) = proot

        ! Create the Cartesian communicator
        call MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, pbc, reorder, cart_comm, ierr)

        ! Get rank within the new communicator
        call MPI_Comm_rank(cart_comm, my_rank, ierr)

        ! Get and store the coordinates of the current MPI task
        call MPI_Cart_coords(cart_comm,my_rank,ndims,my_rank_coords,ierr)

        ! Rank of neighbouring tasks in the four directions
        call MPI_Cart_shift(cart_comm,0,1,my_rank_neighbours(1),my_rank_neighbours(2),ierr)
        
        call MPI_Cart_shift(cart_comm,1,1,my_rank_neighbours(3),my_rank_neighbours(4),ierr)

    end subroutine comms_processor_map


    subroutine comms_finalise()

        integer :: ierr

        call MPI_Finalize(ierr)

    end subroutine comms_finalise

end module comms