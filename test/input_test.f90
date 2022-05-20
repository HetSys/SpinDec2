module input_tester
    use input_params
    use checkpointing
    use netcdf
    use iso_fortran_env

    implicit none

contains

    subroutine input_test()

        integer:: nx, ny, cint, random_seed, err, i, count, ncerr, use_input, current_time
        character(len = 128):: cpi, cpo, problem
        real(kind = real64):: c_min, c_max, m1, m2, k, bfe, t, delta_t, df_tol, ea, eb, tmax, tmin, stab
        real(kind = real64), dimension(:), allocatable:: coeffs
        real(kind = real64), dimension(:, :, :), allocatable:: datas
        real(kind = real64), dimension(:, :), allocatable:: mu
        real(kind = real64), dimension(:), allocatable:: ftot

        call read_params("../input_test.txt",problem, c_min, c_max, coeffs, nx, &
                         ny, m1, m2, ea, eb, Tmin, Tmax, k, bfe, cint, cpi, cpo, &
                          t, delta_t, df_tol, stab, random_seed, use_input, err)

        if (err == -1) then
            print *, "Input test 1 Failed"
            stop
        end if

        if (abs(c_max-0.7) < c_max/1e7 .and. &
            abs(c_min-0.1) < c_min/1e7 .and. &
            nx == 50 .and. &
            ny == 51 .and. &
            abs(m1-1) < m1/1e7 .and. &
            abs(m2-12) < m2/1e7 .and. &
            abs(k-1.6) < k/1e7 .and. &
            abs(bfe-1.7) < bfe/1e7 .and. &
            cint == 50000 .and. &
            cpo == "Out.cp" .and. &
            size(coeffs) == 4 .and. &
            abs(coeffs(1) - 2.1) < coeffs(1) / 1e7 .and. &
            abs(coeffs(2) - 2.2) < coeffs(2) / 1e7 .and. &
            abs(coeffs(3) - 2.3) < coeffs(3) / 1e7 .and. &
            abs(coeffs(4) - 2.4) < coeffs(4) / 1e7 .and. &
            abs(t-1e-7) < t/1e7 .and. &
            abs(delta_t-1e-12) < delta_t/1e7 .and. &
            abs(df_tol-2) < df_tol/1e7 .and. &
            Random_seed == 12345356 .and. &
            problem == "Temp" .and. &
            abs(ea-0.1) < 1/1e7 .and. &
            abs(eb-0.2) < 1/1e7 .and. &
            abs(tmin  - 800) < abs(tmin/1e7) .and. &
            abs(tmax  - 1000) < abs(tmax/1e7) .and. &
            abs(stab-1) < abs(stab/1e7) .and. &
            Use_input == 0) then
        else
            print *, "Input test 1 Failed"
            stop
        end if

        call read_params("../input_test2.txt",problem, c_min, c_max, coeffs, nx, &
                         ny, m1, m2, ea, eb, Tmin, Tmax, k, bfe, cint, cpi, cpo, &
                          t, delta_t, df_tol, stab, random_seed, use_input, err)

        if (err == -1) then
            print *, "Input test 2 Failed"
            stop
        end if

        if (abs(c_max-0.7) < c_max/1e7 .and. &
            abs(c_min-0.1) < c_min/1e7 .and. &
            nx == 50 .and. &
            ny == 51 .and. &
            abs(m1-1) < m1/1e7 .and. &
            abs(m2-12) < m2/1e7 .and. &
            abs(k-1.6) < k/1e7 .and. &
            abs(bfe-1.7) < bfe/1e7 .and. &
            cint == 50000 .and. &
            cpo == "Checkpoint.cpf" .and. &
            cpi == "" .and. &
            size(coeffs) == 5 .and. &
            abs(coeffs(1) - 0) < 1/1e7 .and. &
            abs(coeffs(2) - 0) < 1/1e7 .and. &
            abs(coeffs(3) - 1.7) < abs(coeffs(3) / 1e7) .and. &
            abs(coeffs(4) + 3.4) < abs(coeffs(4) / 1e7) .and. &
            abs(coeffs(5) - 1.7) < abs(coeffs(5) / 1e7) .and. &
            abs(t+1) < abs(t/1e7) .and. &
            abs(delta_t+1) < abs(delta_t/1e7) .and. &
            abs(df_tol+1) < abs(df_tol/1e7) .and. &
            Random_seed == -1 .and. &
            Use_input == 0) then

        else
            print *, "Input test 2 Failed"
            stop
        end if

        call read_params("../input_test3.txt",problem, c_min, c_max, coeffs, nx, &
                         ny, m1, m2, ea, eb, Tmin, Tmax, k, bfe, cint, cpi, cpo, &
                          t, delta_t, df_tol, stab, random_seed, use_input, err)

        if (err /= -1) then
            print *, "Input test 3 Failed"
            stop
        end if

        print *, "All input test passed"

    end subroutine input_test

end module input_tester
