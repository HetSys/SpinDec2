module checkpoint_tester
    use input_params
    use checkpointing
    use netcdf
    use iso_fortran_env
    use grid

    implicit none

contains

    subroutine checkpoint_test()

        integer :: nx, ny, cint, random_seed, err, i, count, ncerr, use_input, &
                    current_time, lastw,wf,wfc,singl,wc
        character(len=128) :: cpi, cpo,problem
        real(kind=real64) :: initial_conc,c_max, c_min, m1, m2, k, bfe, t, delta_t, df_tol,ea,eb,tmin,tmax,stab
        real(kind=real64), dimension(:), allocatable :: coeffs
        real(kind=real64), dimension(:, :, :), allocatable :: datas
        real(kind=real64), dimension(:, :), allocatable :: mu,Tg
        real(kind=real64), dimension(:), allocatable :: ftot

        allocate (datas(1, 1, 1))
        allocate (mu(1, 1))
        allocate (Tg(1, 1))
        allocate (ftot(1))
        allocate (local_grid_conc(1,1))

        datas(1, 1, 1) = 1
        mu(1, 1) = 1
        Tg(1,1) = 1
        ftot(1) = 1
        current_time = 0
        count = 0
        initial_conc = (c_min+c_max)/2
        lastw = 1
        wfc =1
        wc =1

        open (unit=13, file="Out.cp")
        close (13)
        call read_params("../checkpoint_test.txt",problem, c_min, c_max, coeffs, nx, &
                         ny, m1, m2,ea,eb,Tmin,Tmax, k, bfe, cint, cpi, cpo,&
                          t, delta_t, df_tol,stab, random_seed, use_input, err,singl,wf)

        do while (current_time < t / delta_t)
            if (count >= cint) then
                call write_checkpoint_file(datas, mu,Tg, ftot,problem, coeffs, cpo, initial_conc &
                                           , nx, ny, m1, m2, k, bfe, Cint, t, delta_t, current_time,&
                                            df_tol, random_seed,lastw,wf,wfc,singl,wc, ncerr)
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
            Tg(1,1) = Tg(1,1)+4
            if (current_time == 50023) then
                exit
            end if

        end do

         if (cpi /= "") then
             call read_checkpoint_metadata(cpi, problem,initial_conc, coeffs, nx, &
                                     ny, m1, m2, k, bfe, cint, cpo, t, delta_t, df_tol, current_time, &
                                     random_seed, use_input,lastw,wf,wfc,singl,wc, ncerr)
             if (ncerr /= nf90_noerr) then
                 print *, "There was an error reading the checkpoint metadata."
                 print *, "Checkpoint test failed"
                 stop
             end if

             call read_checkpoint_data(datas, mu,Tg, ftot, cpi, use_input, ncerr)
             if (ncerr /= nf90_noerr) then
                 print *, "There was an error reading the checkpoint data."
                 print *, "Checkpoint test failed"
                 stop
             end if
         end if

        if (current_time == 50000 .and. &
            abs(mu(1, 1) - 100001) < mu(1, 1) / 1e7 .and. &
            abs(ftot(1) - 150001) < ftot(1) / 1e7 .and. &
            abs(local_grid_conc(1,1) - 50001) < local_grid_conc(1,1) / 1e7 .and. &
            abs(Tg(1, 1) - 200001) < Tg(1, 1) / 1e7) then
        else
            print *, "Checkpoint test failed"
            stop
        end if

        print *, "All checkpoint tests passed"
    end subroutine checkpoint_test
end module checkpoint_tester
