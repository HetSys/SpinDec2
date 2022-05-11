module checkpoint_tester
    use input_params
    use checkpointing
    use netcdf
    use iso_fortran_env

    implicit none

contains

    subroutine checkpoint_test()

        integer :: nx, ny, cint, random_seed, err, i, count, ncerr, use_input, current_time
        character(len=128) :: cpi, cpo
        real(kind=real64) :: initial_conc, conc_std, m1, m2, k, bfe, t, delta_t, df_tol
        real(kind=real64), dimension(:), allocatable :: coeffs
        real(kind=real64), dimension(:, :, :), allocatable :: datas
        real(kind=real64), dimension(:, :), allocatable :: mu
        real(kind=real64), dimension(:), allocatable :: ftot

        allocate (datas(1, 1, 1))
        allocate (mu(1, 1))
        allocate (ftot(1))

        datas(1, 1, 1) = 1
        mu(1, 1) = 1
        ftot(1) = 1
        current_time = 0
        count = 0

        open (unit=13, file="Out.cp")
        close (13)

        call read_params("../checkpoint_test.txt", initial_conc, conc_std, coeffs, nx, &
                         ny, m1, m2, k, bfe, cint, cpi, cpo, t, delta_t, df_tol, random_seed, use_input, err)

        do while (current_time < t / delta_t)
            if (count >= cint) then
                call write_checkpoint_file(datas, mu, ftot, coeffs, cpo, initial_conc, conc_std &
                                           , nx, ny, m1, m2, k, bfe, Cint, t, delta_t, current_time, df_tol, random_seed, ncerr)
                if (ncerr /= nf90_noerr) then
                    print *, "There was an error writing the checkpoint file."
                    print *, "Checkpoint test failed"
                    stop
                end if
                count = 0
            end if
            count = count + 1
            current_time = current_time + 1
            datas(1, 1, 1) = datas(1, 1, 1) + 1
            mu(1, 1) = mu(1, 1) + 2
            ftot(1) = ftot(1) + 3
            if (current_time == 50023) then
                exit
            end if

        end do

        if (cpi /= "") then
            call read_checkpoint_in(datas, mu, ftot, cpi, initial_conc, conc_std, coeffs, nx, &
                                    ny, m1, m2, k, bfe, cint, cpo, t, delta_t, df_tol, current_time, random_seed, use_input, ncerr)
            if (ncerr /= nf90_noerr) then
                print *, "There was an error reading the checkpoint file."
                print *, "Checkpoint test failed"
                stop
            end if
        end if

        if (current_time == 50000 .and. &
            abs(mu(1, 1) - 100001) < mu(1, 1) / 1e7 .and. &
            abs(ftot(1) - 150001) < ftot(1) / 1e7 .and. &
            abs(datas(1, 1, 1) - 50001) < datas(1, 1, 1) / 1e7) then
        else
            print *, "Checkpoint test failed"
            stop
        end if

        print *, "All checkpoint tests passed"
    end subroutine checkpoint_test
end module checkpoint_tester
