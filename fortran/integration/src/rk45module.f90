module rk45module

    use kindmodule

    implicit none

    ! Input/Output for integrator 
    type IntegratorInOut
        integer :: n ! State size
        complex(wp) :: t0 ! Initial time
        complex(wp) :: tf ! Final time
        complex(wp) :: h0=1.0D-4 ! Initial step size
        integer :: stepmax = 10000  ! Maximum number of steps
        complex(wp) :: hmin = 1.0D-6 ! Minimum step size
        real(wp) :: tol = 1.0D-10 ! Tolerance
        complex(wp), dimension(:), allocatable :: x0 ! Initial state 
        complex(wp), dimension(:), allocatable :: tout ! Time output 
        complex(wp), dimension(:), allocatable :: eout ! Error output 
        complex(wp), dimension(:), allocatable :: hout ! Step size output 
        complex(wp), dimension(:,:), allocatable :: xout ! State output 
    end type IntegratorInOut

    ! Runge-Kutta-Fehlberg coefficients
    complex(wp), parameter :: alpha(5) = &
        [1.0/4.0, 3.0/8.0, 12.0/13.0, 1.0, 1.0/2.0]

    complex(wp), parameter :: beta(5,5) = transpose(reshape([ &
        1.0/4.0, 0.0, 0.0, 0.0, 0.0, &
        3.0/32.0, 9.0/32.0, 0.0, 0.0, 0.0, &
        1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0, 0.0, 0.0, &
        439.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0, 0.0, &
        -8.0/27.0, 2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0], &
        shape(beta), order=[2,1]))

    complex(wp), parameter :: c(5) = &
        [25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0]

    complex(wp), parameter :: chat(6) = &
        [16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, &
        -9.0/50.0, 2.0/55.0]

contains 

    subroutine rk45(f_ptr, io)

        use kindmodule
        use derivkepmodule

        ! DECLARATION
        implicit none

        ! In/Out
        type(IntegratorInOut), intent(inout) :: io 
        
        abstract interface
        function func(t, x)
        use kindmodule
        complex(wp),intent(in) :: t
        complex(wp),dimension(6),intent(in) :: x
        complex(wp),dimension(6) :: func 
        end function func
        end interface 
        
        procedure (func), pointer :: f_ptr => null ()  ! Derivative function 

        ! Local
        complex(wp) :: t ! Time
        complex(wp),dimension(:,:),allocatable :: f ! Derivatives
        complex(wp) :: h ! Step size
        complex(wp),dimension(:),allocatable :: x0
        complex(wp),dimension(:),allocatable :: xhatadd
        complex(wp),dimension(:),allocatable :: xadd
        complex(wp),dimension(:),allocatable :: xhat
        complex(wp),dimension(:),allocatable :: x
        complex(wp),dimension(:,:),allocatable :: xstore
        complex(wp),dimension(:),allocatable :: tstore
        complex(wp),dimension(:),allocatable :: estore
        complex(wp),dimension(:),allocatable :: hstore
        integer :: step
        complex(wp) :: error
        integer :: i

        allocate(f(6,io%n))
        allocate(x0(io%n))
        allocate(xhatadd(io%n))
        allocate(xadd(io%n))
        allocate(xhat(io%n))
        allocate(x(io%n))
        allocate(xstore(io%stepmax,io%n))
        allocate(tstore(io%stepmax))
        allocate(estore(io%stepmax))
        allocate(hstore(io%stepmax))

        ! EXECUTION

        ! Initialization
        x = io%x0
        h = io%h0
        step = 1

        do
            if ((real(t) >= real(io%tf)) .OR. (step > io%stepmax)) then
                exit
            end if

            ! Evaluate derivatives
            f(1,:) = h*f_ptr(t, x)
            f(2,:) = h*f_ptr(t + alpha(1)*h, x + beta(1,1)*f(1,:))
            f(3,:) = h*f_ptr(t + alpha(2)*h, x + beta(2,1)*f(1,:) &
                + beta(2,2)*f(2,:))
            f(4,:) = h*f_ptr(t + alpha(3)*h, x + beta(3,1)*f(1,:) &
                + beta(3,2)*f(2,:) + beta(3,3)*f(3,:))
            f(5,:) = h*f_ptr(t + alpha(4)*h, x + beta(4,1)*f(1,:) & 
                + beta(4,2)*f(2,:) + beta(4,3)*f(3,:) + beta(4,4)*f(4,:))
            f(6,:) = h*f_ptr(t + alpha(5)*h, x + beta(5,1)*f(1,:) &
                + beta(5,2)*f(2,:) + beta(5,3)*f(3,:) + beta(5,4)*f(4,:) &
                + beta(5,5)*f(5,:))

            ! Evaluate integrated states
            xhatadd = 0
            xadd = 0
            do i=1,6
                xhatadd = xhatadd + chat(i)*f(i,:)
            end do
            do i=1,5
                xadd = xadd + c(i)*f(i,:)
            end do
            xhat = x + xhatadd
            x = x + xadd

            ! Check difference between high/low order 
            error = norm2(real(x - xhat))

            ! Save current step
            xstore(step,:) = x
            tstore(step) = t
            estore(step) = error
            hstore(step) = h

            ! Step forward (TODO step size control)
            t = t + h
            step = step + 1

        end do 

        ! Populate output
        allocate(io%xout(step-1,io%n))
        allocate(io%tout(step-1))
        allocate(io%eout(step-1))
        allocate(io%hout(step-1))

        do i=1,step-1
            io%xout(i,:) = xstore(i,:)
            io%tout(i) = tstore(i) 
            io%eout(i) = estore(i)
            io%hout(i) = hstore(i) 
        end do

        ! Cleanup
        deallocate(f)
        deallocate(x0)
        deallocate(xhatadd)
        deallocate(xadd)
        deallocate(xhat)
        deallocate(x)
        deallocate(xstore)
        deallocate(tstore)
        deallocate(estore)
        deallocate(hstore)

    end subroutine rk45

end module rk45module
