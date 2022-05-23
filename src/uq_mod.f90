module uq

    use iso_fortran_env
    use cahn_hilliard
    use grid
    use potentials
    use spectral
    use io
    use free_energy
    use input_params
    use checkpointing

    implicit none
contains
    subroutine run_sim(fn,exp,dt_in,Nt_max)
        real(real64), dimension(:, :, :), allocatable :: c ! conc. grid
        real(real64), dimension(:, :), allocatable :: c_new, c_out ! new conc. grid
        real(real64), dimension(:, :), allocatable :: mu   ! bulk chem. pot.
        real(real64), dimension(:, :), allocatable :: Q    ! total chem. pot.
        real(real64), dimension(:, :), allocatable :: dQ   ! 2nd derivative of Q
        real(real64), dimension(:, :), allocatable :: M    ! Mobility field
        real(real64), dimension(:, :), allocatable :: T ! Temp
        real(real64), dimension(:), allocatable :: a ! user inputted polynomial coefficients
        real(real64), dimension(:), allocatable :: F_tot ! Total free energy with time
        real(real64), dimension(:, :), allocatable :: f_b
        real(real64) :: c0, c_std !initial grid mean and std sample - if using normal dist
        real(real64) :: c_min, c_max !initial grid lower and upper boubds - for default uniform dist
        real(real64) :: T_min, T_max !initial temp grid mean and std sample
        real(real64) :: dx, dy, dt ! spatial and temporal grid spacings
        real(real64) :: Kappa ! free energy gradient parameter
        real(real64) :: t_end !end time
        real(real64) :: MA, MB ! Mobility's
        real(real64) :: EA, EB ! exciation energy
        real(real64) :: stab !Stabilization_Term
        real(real64) :: bfe, df_tol!Placeholder (These were in the input file but df_tol hasn't been used in any code)
        integer :: Nx, Ny, Nt, Nc
        integer :: k, count,thread ! counters
        integer :: cint, random_seed, err, use_input, current_iter, ncerr !checkpointing_interval, random seed,error var
        character(len=128) :: cpi, cpo ! checkpointing files
        character(len=128) :: problem
        character(len=*),intent(in) :: fn
        logical,intent(out) :: exp
        real(real64) ,intent(in) :: dt_in
        integer, intent(in) :: Nt_max



        ! Only run files in test for now
        call read_params(fn,problem, c_min, c_max, a, nx, &
                         ny, ma, mb,ea,eb,T_min,T_max, kappa, bfe, cint, cpi, &
                         cpo, t_end, dt, df_tol,stab, random_seed, use_input, err)

        if (err == -1) then
            print *, "There was an issue with the input file please check and try again"
            stop
        end if

        current_iter = 2
        ! come back to this

        c0 = (c_min+c_max)/2
        c_std = 0

        ! Allocate T grid
        if (.not. allocated(T)) then
            allocate (T(Nx, Ny))
            T = 0.0
        end if

        if (problem == 'Temp') then
            call grid_init(T, Nx, Ny, T_min, T_max)
        end if


        ! Set seed
        call get_seed(random_seed)
        !print*, random_seed
        dt = dt_in
        Nt = Nt_max
        dx = 1.0_real64/real(nx)
        dy = 1.0_real64/real(ny)

        if( problem /= 'Spectral') then
            !if (dt > min(0.1*dx**4, 0.1*dy**4)) then
            !    print*, 'Warning time-step unstable, setting to default stable value'
            !    dt = min(0.1*dx**4, 0.1*dy**4)
            !end if
        end if


        !Find polynomial coefficients size
        Nc = size(a)
        ! Allocate c grid
        if (.not. allocated(c)) then
            allocate (c(Nx, Ny, 3))
            c = 0.0
        end if

         ! Allocate M grid
        if (.not. allocated(M)) then
            allocate (M(Nx, Ny))
            M = 0.0
        end if

        ! Allocate mu grid
        if (.not. allocated(mu)) then
            allocate (mu(Nx, Ny))
            mu = 0.0
        end if

        ! Allocate F_tot
        if (.not. allocated(F_tot)) then
            allocate (F_tot(Nt))
            F_tot = 0.0
        end if

         ! Allocate Q grid
        if (.not. allocated(Q)) then
            allocate (Q(Nx, Ny))
            Q = 0.0
        end if

         ! Allocate c_new grid
        if (.not. allocated(c_new)) then
            allocate (c_new(Nx, Ny))
            c_new = 0.0
        end if

        if (.not. allocated(c_out)) then
            allocate (c_out(Nx, Ny))
            c_out = 0.0
        end if
        ! Initialize c grid

        call grid_init(c(:, :, 2), Nx, Ny, c_min, c_max)

        if (problem == 'Temp') then
            call grid_init(T, Nx, Ny, T_min, T_max)
        end if

        ! Get Initial Bulk Free Energy over space
        call bulk_free_energy(f_b, c(:, :, 2), a)

        ! Calculate Initial F(t)
        call total_free_energy(F_tot(1), c(:, :, 2), f_b, dx, dy, kappa)


        deallocate (f_b)

        ! Grid evolution
        do k = current_iter, Nt


            ! Get bulk chemical potentials
            call bulk_potential(mu, c(:, :,2), a)

            ! Get total chemical potentials
            call total_potential(Q, mu, c(:, :, 2), dx, dy, Kappa)

            ! Get Mobility Field
            call Mobility(M,MA,MB, EA, EB, c0, c(:, :, 2), T, problem)

            !print*, M(1,1), M(1,2), M(6,7)

            ! Get new concentrations for current timesteps

            if (problem == 'Spectral') then
                if(k == 2) then
                    call spectral_method_iter(c(:, :, 2),c(:, :, 2),a,dt,M(1,1),Kappa,c_out,1,stab)
                    c(:,:,3) = c_out
                else
                    call spectral_method_iter(c(:, :, 2),c(:, :, 1),a,dt,M(1,1),Kappa,c_out,0,stab)
                    c(:,:,3) = c_out
                end if
            else
                call time_evoloution_new(c(:, :, 2),c_new,M,Q,dx,dy,dt, Nx, Ny)

                c(:, :, 3) = c_new(:, :)
            end if




            ! Get Bulk Free Energy over space
            call bulk_free_energy(f_b, c(:, :, 3), a)

            ! Calculate F(t)
            call total_free_energy(F_tot(k), c(:, :, 3), f_b, dx, dy, kappa)

            deallocate (f_b)

            c(:,:,1) = c(:,:,2)
            c(:,:,2) = c(:,:,3)

            !print*, maxval(c(:,:,3))

            if(maxval(c(:,:,3)) > 1.1 .or. minval(c(:,:,3))<-0.1) then
                exp = .true.
                return
            else
                exp = .false.
            end if


        end do

    end subroutine run_sim

    subroutine find_critical_timestep(fn,dt_ci,dt_si,Nt_max,level,ts)
        real(real64) ,intent(in) :: dt_ci,dt_si
        integer,intent(in) :: Nt_max,level
        character(len=*),intent(in) :: fn
        integer :: i,l
        logical :: exp
        real(real64),intent(out) :: ts
        real(real64) ::dt_c,dt_s
        dt_c=dt_ci
        dt_s=dt_si
        do l=1,level
            i = 0
            exp = .false.
            do while(exp .eqv. .false.)
                call run_sim(fn,exp,dt_s+dt_c*i,Nt_max)
                i = i + 1
                if(dt_s+dt_c*i > 10) then
                    exp = .true.
                end if
            end do

            dt_s = dt_s+dt_c*(i-2)
            dt_c = dt_c/10.0_real64
            !print*,dt_s

        end do

        ts=dt_s
        print*, dt_s
    end subroutine find_critical_timestep


end module uq
