!Temp_max < temp_min
module input_params

    use iso_fortran_env

    implicit none

contains

    function read_poly(var, err) result(coeffs)
        ! A function that reads a string of comma seperated values into an array

        character(*), intent(in) :: var
        integer, intent(inout) :: err
        character(len=128) :: a, b
        integer :: lens, i, j, read_counter, ierr
        real(kind=real64), dimension(:), allocatable :: coeffs
        read_counter = 0
        err = 0

        do i = 1, len(var)
            read_counter = read_counter + 1
            if (var(i:i) == ",") then
                exit
            end if
        end do

        b = trim(adjustl(var(read_counter + 1:)))
        a = trim(adjustl(var(:read_counter - 1)))
        read (var, *, iostat=ierr) lens

        if (ierr /= 0) then
            print *, "An error occured reading f(c). &
            &Please check input file and try again"
            err = -1
            return
        end if

        if (lens < 0) then
            err = -1
            print *, "You must have at least 1 coefficent. Please check your &
            &your input file and try again"
            allocate (coeffs(1))
            return
        end if

        allocate (coeffs(lens))
        ! print*,b

        do i = 1, lens
            read_counter = 0
            do j = 1, len(b)
                read_counter = read_counter + 1
                if (b(j:j) == ",") then
                    exit
                end if
            end do

            a = trim(adjustl(b(:read_counter - 1)))
            b = trim(adjustl(b(read_counter + 1:)))
            read (a, *, iostat=ierr) coeffs(i)

            if (ierr /= 0) then
                print *, "An error occured reading f(c). &
                &Please check input file and try again"
                err = -1
                return

            end if
        end do

    end function read_poly

    subroutine read_params(fn,prob, conc_min, conc_max, coeffs, Nx, Ny, M1, M2,EA,EB,temp_min,temp_max, k, bfe, &
                           Cint, cpi, cpo, t, delta_t, df_tol,stab, random_seed, use_input, err)

        integer, parameter :: infile = 15
        character(len=128) :: name, var
        character(*), intent(in) :: fn
        integer :: read_counter, i, iostat, ierr
        integer, intent(out) :: nx, ny, cint, random_seed, use_input
        integer, intent(out) :: err
        character(len=128), intent(out) :: cpi, cpo,prob
        real(kind=real64), intent(out) :: conc_min, conc_max, m1, m2, k, bfe, t, &
                                            delta_t, df_tol,ea,eb,temp_max,temp_min,stab
        real(kind=real64), intent(out), dimension(:), allocatable :: coeffs
        logical :: f_exists
        integer :: req, req2, cofs,treq,sreq

        ! Setting deafault values for optional params
        req = 0
        req2 = 0
        treq = 0
        sreq = 0
        cofs = 0
        err = 0
        random_seed = -1
        delta_t = -1
        t = -1
        cpi = ""
        cpo = "checkpoint.cpf"
        use_input = 0
        df_tol = -1
        ea = -1
        eb = -1
        temp_max = -1
        temp_min = -1
        stab = -1

        ! reading file
        ! There is a large amount of if statements below to ensure values are
        ! sensisble (no negative max time etc)
        open (unit=infile, file=fn, status='old', action="read", iostat=ierr)
        if (ierr /= 0) then
            print *, "An error occured when reading the input file&
             & (It may be missing). Please check and try again."
            err = -1
        end if

        do
            if (err == -1) then
                exit
            end if

            read (infile, "(a)", iostat=iostat) name
            if (iostat < 0) then
                exit
            end if

            !print*,name
            read_counter = 0
            do i = 1, len(name)
                read_counter = read_counter + 1
                if (name(i:i) == "=") then
                    exit
                end if
            end do

            var = trim(adjustl(name(read_counter + 1:)))
            name = trim(adjustl(name(:read_counter - 1)))

            if (name == "Concentration_min") then
                req = req + 1
                read (var, *, iostat=ierr) conc_min

                if (ierr /= 0) then
                    print *, "An error occured reading Concentration_min. &
                    &Please check input file and try again"
                    err = -1
                    exit
                end if

            else if (name == "Concentration_max") then
                req = req + 1
                read (var, *, iostat=ierr) conc_max
                if (ierr /= 0) then
                    print *, "An error occured reading Concentration_max. &
                    &Please check input file and try again"
                    err = -1
                    exit
                end if

            else if (name == "f(c)") then
                req2 = req2 + 1
                cofs = cofs + 1
                coeffs = read_poly(var, err)

                if (err == -1) then
                    deallocate (coeffs)
                end if

            else if (name == "Domain_x_size") then
                req = req + 1
                read (var, *, iostat=ierr) nx

                if (ierr /= 0) then
                    print *, "An error occured reading Domain_x_size. &
                    &Please check input file and try again"
                    err = -1
                    exit
                end if

                if (nx < 1) then
                    print *, "Domain_x_size must be >= 1.&
                    & Please check you input file and try again"
                    err = -1
                end if

            else if (name == "Domain_y_size") then
                req = req + 1
                read (var, *, iostat=ierr) ny

                if (ierr /= 0) then
                    print *, "An error occured reading Domain_y_size. &
                    &Please check input file and try again"
                    err = -1
                    exit
                end if

                if (ny < 1) then
                    print *, "Domain_y_size must be >= 1.&
                    & Please check you input file and try again"
                    err = -1
                end if

            else if (name == "Mobility_A") then
                req = req + 1
                read (var, *, iostat=ierr) m1

                if (ierr /= 0) then
                    print *, "An error occured reading Mobility_A. &
                    &Please check input file and try again"
                    err = -1
                    exit
                end if

                if (m1 < 0) then
                    print *, "Mobility_A must be >= 0.&
                    & Please check you input file and try again"
                    err = -1
                end if

            else if (name == "Mobility_B") then
                req = req + 1
                read (var, *, iostat=ierr) m2

                if (ierr /= 0) then
                    print *, "An error occured reading Mobility_B. &
                    &Please check input file and try again"
                    err = -1
                    exit
                end if

                if (m2 < 0) then
                    print *, "Mobility_B must be >= 0.&
                    & Please check you input file and try again"
                    err = -1
                end if

            else if (name == "free_energy_gradient_parameter") then
                req = req + 1
                read (var, *, iostat=ierr) k

                if (ierr /= 0) then
                    print *, "An error occured reading &
                    &free_energy_gradient_parameter. &
                    &Please check input file and try again"
                    err = -1
                    exit
                end if

            else if (name == "Bulk_free_energy") then
                req2 = req2 + 1
                read (var, *, iostat=ierr) bfe

                if (ierr /= 0) then
                    print *, "An error occured reading Bulk_free_energy. &
                    &Please check input file and try again"
                    err = -1
                    exit
                end if

            else if (name == "Checkpointing_interval") then
                req = req + 1
                read (var, *, iostat=ierr) cint

                if (ierr /= 0) then
                    print *, "An error occured reading Checkpointing_interval. &
                    &Please check input file and try again"
                    err = -1
                    exit
                end if

                if (cint < 1) then
                    print *, "Checkpointing_interval must be >= 1.&
                    & Please check you input file and try again"
                    err = -1
                end if

            else if (name == "Checkpoint_input_file") then
                read (var, *, iostat=ierr) cpi

                if (ierr /= 0) then
                    print *, "An error occured reading Checkpoint_input_file. &
                    &Please check input file and try again"
                    err = -1
                    exit
                end if

                inquire (file=cpi, exist=f_exists)
                if (f_exists .eqv. .false.) then
                    print *, "The checkpoint input file does not exist.&
                    &  Plese check input file and try again"
                    err = -1
                end if

            else if (name == "Checkpoint_output_file") then
                read (var, *, iostat=ierr) cpo

                if (ierr /= 0) then
                    print *, "An error occured reading Checkpoint_output_file. &
                    &Please check input file and try again"
                    err = -1
                    exit
                end if

            else if (name == "Max_time") then
                read (var, *, iostat=ierr) t

                if (ierr /= 0) then
                    print *, "An error occured reading Max_time. &
                    &Please check input file and try again"
                    err = -1
                    exit
                end if

                if (t < 0) then
                    print *, "Max_time must be >= 0.&
                    & Please check you input file and try again"
                    err = -1
                end if

            else if (name == "time_step") then
                read (var, *, iostat=ierr) delta_t

                if (ierr /= 0) then
                    print *, "An error occured reading time_step. &
                    &Please check input file and try again"
                    err = -1
                    exit
                end if

                if (delta_t < 0) then
                    print *, "time_step must be >= 0.&
                    & Please check you input file and try again"
                    err = -1
                end if

            else if (name == "dF_tolerance") then
                read (var, *, iostat=ierr) df_tol

                if (ierr /= 0) then
                    print *, "An error occured reading dF_tolerance. &
                    &Please check input file and try again"
                    err = -1
                    exit
                end if

                if (df_tol < 0) then
                    print *, "dF_tolerance must be >= 0.&
                    & Please check you input file and try again"
                    err = -1
                end if

            else if (name == "Random_seed") then
                read (var, *, iostat=ierr) random_seed

                if (ierr /= 0) then
                    print *, "An error occured reading Random_seed. &
                    &Please check input file and try again"
                    err = -1
                    exit
                end if

                if (random_seed < -1) then
                    print *, "Random_seed must be >= -1.&
                    & Please check you input file and try again"
                    err = -1
                end if

            else if (name == "Use_input") then
                read (var, *, iostat=ierr) use_input

                if (ierr /= 0) then
                    print *, "An error occured reading Random_seed. &
                    &Please check input file and try again"
                    err = -1
                    exit
                end if

                if (use_input /= 0 .and. use_input /= 1) then
                    print *, "use_input must be 1(True) or 0(False).&
                    & Please check you input file and try again"
                    err = -1
                end if
            else if (name == "Exitation_A") then
                treq = treq + 1
                read (var, *, iostat=ierr) ea

                if (ierr /= 0) then
                    print *, "An error occured reading Exitation_A. &
                    &Please check input file and try again"
                    err = -1
                    exit
                end if

                if (ea < 0) then
                    print *, "Exitation_A must be >= 0.&
                    & Please check you input file and try again"
                    err = -1
                end if

            else if (name == "Exitation_B") then
                treq = treq + 1
                read (var, *, iostat=ierr) eb

                if (ierr /= 0) then
                    print *, "An error occured reading Exitation_B. &
                    &Please check input file and try again"
                    err = -1
                    exit
                end if

                if (eb < 0) then
                    print *, "Exitation_B must be >= 0.&
                    & Please check you input file and try again"
                    err = -1
                end if
            else if (name == "Temperature_max") then
                treq = treq + 1
                read (var, *, iostat=ierr) temp_max

                if (ierr /= 0) then
                    print *, "An error occured reading Temperature_max. &
                    &Please check input file and try again"
                    err = -1
                    exit
                end if

                if (temp_max < 0) then
                    print *, "Temperature_max must be >= 0.&
                    & Please check you input file and try again"
                    err = -1
                end if

            else if (name == "Temperature_min") then
                treq = treq + 1
                read (var, *, iostat=ierr) temp_min

                if (ierr /= 0) then
                    print *, "An error occured reading Temperature_min. &
                    &Please check input file and try again"
                    err = -1
                    exit
                end if

                if (temp_min< 0) then
                    print *, "Temperature_min must be >= 0.&
                    & Please check you input file and try again"
                    err = -1
                end if
            else if (name == "Problem") then
                req = req +1
                read (var, *, iostat=ierr) prob

                if (ierr /= 0) then
                    print *, "An error occured reading Checkpoint_input_file. &
                    &Please check input file and try again"
                    err = -1
                    exit
                end if

                if(prob /= "Spectral" .and. prob /= "Temp" .and. prob /= "Constant" .and. prob /= "NonTemp") then
                    print*, "Problem must be one of; Spectral, Temp, Constant, NonTemp"
                    err = -1
                end if
            else if (name == "Stabilization_Term") then
                sreq = sreq +1
                read (var, *, iostat=ierr) stab

                if (ierr /= 0) then
                    print *, "An error occured reading Stabilization_Term. &
                    &Please check input file and try again"
                    err = -1
                    exit
                end if

                if (eb < 0) then
                    print *, "Stabilization_Term must be >= 0.&
                    & Please check you input file and try again"
                    err = -1
                end if

            end if
            !print*, name,var
        end do

        if (req < 9) then
            err = -1
            print *, "Not all required inputs were found. Plase check an try again."
        else if (req2 == 0) then
            err = -1
            print *, "Not all required inputs were found. Plase check an try again."
        else
            if(conc_max < conc_min) then
                err = -1
                print*, "Concentration_max must be >= Concentration_min."
            end if
        end if

        if (req2 == 2) then
            print *, "The bulk free energy will be ignored as there is a custom f(c)"
        end if

        if (cofs == 0) then
            allocate (coeffs(5))
            coeffs(1) = 0
            coeffs(2) = 0
            coeffs(3) = bfe
            coeffs(4) = -2*bfe
            coeffs(5) = bfe
        end if

        if(prob == "Temp") then
            if(treq < 4) then
                err  = -1
                print*, "Not all required inputs found for problem Temp"
            else
                if(temp_max < temp_min) then
                    err = -1
                    print*, "Temperature_max must be >= Temperature_min."
                end if
            end if
        end if

        if(prob == "Spectral") then
            if(sreq < 1) then
                err = -1
                print*, "Not all required inputs found for problem Spectral"
            end if
        end if


    close(infile)
    end subroutine read_params

end module input_params
