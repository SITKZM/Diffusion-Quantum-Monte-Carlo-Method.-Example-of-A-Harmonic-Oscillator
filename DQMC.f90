! gfortran -o DQMC.out Ziggurat.f90 DQMC.f90 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
program DQMC
    use Ziggurat
    implicit none
    ! parameters
    integer, parameter :: N_ini = 500
    integer, parameter :: N_max = 2000
    integer, parameter :: seed_ini = 43
    integer, parameter :: n_b = 101
    integer, parameter :: tau_0 = 1000
    double precision, parameter :: Deltatau = 0.05
    double precision, parameter :: x_min = -5.
    double precision, parameter :: x_max = 5.
    double precision, parameter :: x_0 = 0.
    
    ! variables
    integer :: N_0 = N_ini
    integer :: N_1
    integer :: seed = seed_ini
    double precision :: psips(N_max) = 0.
    double precision :: E_R = 0.1
    double precision :: E_Rs_first(tau_0 + 1) = 0
    double precision :: E_Rs_last(tau_0) = 0
    double precision :: phi_0(n_b) = 0
    double precision :: std_E_R

    ! loop index etc
    integer :: loop, unit_write_result
    logical :: flg = .false.

    call Initialize_replicas()
    E_Rs_first(1) = E_R
    
    do loop = 1, tau_0
        call walk()
        call branch()
        E_Rs_first(loop + 1) = E_R
        if ( N_0 == 0 ) then
            print *, "All particles are annihilated."
            exit
        else if (N_0 >= N_max) then
            print *, "Number of particles is diverged."
            exit
        end if
    end do

    do loop = 1, tau_0
        call walk()
        call branch()
        call count()
        E_Rs_last(loop) = E_R
    end do

    E_R = sum(E_Rs_last) / tau_0
    std_E_R = std(E_Rs_last, tau_0, E_R)
    call write_files()
contains
    subroutine Initialize_replicas()
        integer :: i

        do i = 1, N_0
            psips(i) = x_0
        end do
    end subroutine Initialize_replicas

    subroutine walk()
        integer :: i
        double precision :: rho
        integer :: kn(128)
        real :: fn(128)
        real :: wn(128)

        do i = 1, N_0
            seed = shr3(seed)
            call r4_nor_setup(kn, fn, wn)
            rho = r4_nor(seed, kn, fn, wn)

            psips(i) = psips(i) + rho * sqrt(Deltatau)
        end do
    end subroutine walk

    subroutine branch()
        integer :: i, j
        integer :: sum_m_n
        integer :: m_n
        double precision :: W
        double precision :: u
        double precision :: V
        double precision :: V_average
        double precision :: psips_new(N_max)

        V_average = 0
        N_1 = N_0
        sum_m_n = 0
        psips_new = 0

        ! calculate m_n, <V>_n
        do i = 1, N_0
            V = psips(i)**2 / 2
            W = exp(-(V - E_R) * Deltatau)

            u = r4_uni(seed)
            m_n = min(int(W + u), 3)

            N_1 = N_1 + (m_n - 1)
            if ( N_1 > N_max ) then
                exit
            end if

            if ( m_n > 0 ) then
                do j = sum_m_n + 1, sum_m_n + m_n
                    psips_new(j) = psips(i)
                end do
            end if
            sum_m_n = sum_m_n + m_n

            V_average = V_average + V / N_0
        end do

        if ( N_0 == N_1 ) then
            flg = .true.
        end if

        ! update
        if ( flg ) then
            E_R = V_average
        else
            E_R = V_average + (1 - N_1 / N_0) / Deltatau
        end if
        
        psips = psips_new
        print *, loop, ".", N_0, N_1, E_R, flg
        N_0 = N_1
    end subroutine branch

    subroutine count()
        integer :: i, j
        double precision :: histogram(n_b)
        double precision :: x
        double precision :: box_width = (x_max - x_min) / n_b

        histogram = 0
        do i = 1, N_0
            x = psips(i)
            do j = 1, n_b
                if ( x_min + box_width * (j - 1) < x .and. x < x_min + box_width * j ) then
                    histogram(j) = histogram(j) + 1
                end if
            end do
        end do

        histogram = histogram / sqrt(box_width * sum(histogram**2))
        phi_0 = phi_0 + histogram / tau_0
    end subroutine count

    function std(arr, dim, average) result(result)
        integer :: i
        integer :: dim
        double precision :: arr(dim)
        double precision :: average
        double precision :: result

        result = 0

        do i = 1, dim
            result = result + (arr(i) - average)**2 / dim
        end do

        result = sqrt(result)
    end function std

    subroutine write_files()
        integer :: i
        double precision :: box_width = (x_max - x_min) / n_b
        double precision :: ledge, redge

        open(newunit = unit_write_result, file = "info.txt")
            write(unit_write_result, *) "Input"
            write(unit_write_result, *) "N_0 =", N_ini
            write(unit_write_result, *) "N_max =", N_max
            write(unit_write_result, *) "seed =", seed_ini
            write(unit_write_result, *) "tau_0 =", tau_0
            write(unit_write_result, *) "Deltatau =", Deltatau
            write(unit_write_result, *) "xmin, xmax =", x_min, x_max
            write(unit_write_result, *) "n_b =", n_b
            write(unit_write_result, *) "x_0 =", x_0
            write(unit_write_result, *) "Output"
            write(unit_write_result, *) "E_0 =", E_R
            write(unit_write_result, *) "std_E_R =", std_E_R
            write(unit_write_result, *) "time evolution of E_R: ER_evolution.txt"
            write(unit_write_result, *) "section enclosed by the box, phi_0(x): phi_0.txt"
        close(unit_write_result)

        open (newunit=unit_write_result, file="ER_evolution.txt")
            do i = 1, tau_0 + 1
                write (unit_write_result, *) Deltatau * (i - 1), E_Rs_first(i)
            end do
        close (unit_write_result)

        open (newunit=unit_write_result, file="phi_0.txt")
            do i = 1, n_b
                ledge = x_min + box_width * (i - 1)
                redge = x_min + box_width * i
                write (unit_write_result, *) ledge, redge, phi_0(i)
            end do
        close (unit_write_result)
    end subroutine write_files
end program DQMC