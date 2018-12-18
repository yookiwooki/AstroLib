module rk45

    use kindmodule

    implicit none

    ! Input/Output for integrator 
    type IntegratorInOut
        integer :: n ! State size
        integer :: stepmax = 100000 ! Maximum number of steps
        external :: func ! Derivative function 
        complex(wp) :: t0 ! Initial time
        real(wp) :: tol = 1.0D-10 ! Tolerance
        complex(wp) :: dtmin = 1.0D-4 ! Minimum timestep size
        complex(wp) :: tf ! Final time
        complex(wp) :: x0(n) ! Initial state 
        complex(wp) :: tout(stepmax) ! Integrator time output 
        complex(wp), allocatable :: xout(stepmax,n) ! Integ. state out
    end type IntegratorInOut

    ! Runge-Kutta-Fehlberg coefficients
    complex(wp), parameter :: alpha(5) = &
        [1.0/4.0, 3.0/8.0, 12.0/13.0, 1.0, 1.0/2.0]

    complex(wp), parameter :: beta(5,5) = transpose(reshape([ &
        1.0/4.0, 0.0, 0.0, 0.0, 0.0, &
        3.0/32.0, 9.0/32.0, 0.0, 0.0, 0.0, &
        1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0, 0.0, 0.0, &
        439.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0, 0.0, &
        -8.0/27.0, 2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0],
    shape(beta)))

    complex(wp), parameter :: c(5) = &
        [25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0]

    complex(wp), parameter :: chat(6) = &
        [16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, &
        -9.0/50.0, 2.0/55.0]

contains subroutine rk45(io)

    use kindmodule

    ! DECLARATION
    implicit none
    ! In/Out
    type(IntegratorInOut), intent(inout) :: io 
    ! Local

    ! EXECUTION

end subroutine rk45

end module rk45
