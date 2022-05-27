program main

    use iso_fortran_env
    use uq

    real(real64) :: dt_c,dt_s,ts
    integer :: i,number
    character(len=1024) :: filename

    number = 14

    thread =  fftw_init_threads()

    dt_c = 1e-2_real64
    dt_s = 0_real64

    open (unit = number, file = "output.txt")

    do i=0,9999
        write (filename, "(A6,I0,A4)") "input_", i,".txt"
        print *, trim(filename)
        call find_critical_timestep(trim(filename),dt_c,dt_s,10,7,ts)
        write(number,*) ts
    end do

    close(number)



end program main
