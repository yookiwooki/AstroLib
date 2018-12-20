! rk45module.f90 - Runge-Kutta-Fehlberg ordinary differential equation solver
module rk45module

    use kindmodule

    implicit none

    ! Input/Output for integrator 
    type IntegratorIn
        integer :: n                                     ! State size
        complex(wp) :: t0                                ! Initial time
        complex(wp) :: tf                                ! Final time
        complex(wp) :: h0=1.0D-5                         ! Initial step size
        integer :: stepmax = 100000                      ! Max number of steps
        real(wp) :: hmin = 1.0D-6                     ! Min step size
        real(wp) :: tol = 1.0D-10                        ! Tolerance
        complex(wp), dimension(:), allocatable :: x0     ! Initial state 
        real(wp) :: rparam = 0.5  ! Step size reduce tuning
        real(wp) :: iparam = 1.5  ! Step size increase tuning
        real(wp) :: dbparam = 0.9 ! Deadband for increasing step size
    end type IntegratorIn

    type IntegratorOut
        complex(wp), dimension(:), allocatable :: tout   ! Time output 
        complex(wp), dimension(:), allocatable :: eout   ! Error output 
        complex(wp), dimension(:), allocatable :: hout   ! Step size output 
        complex(wp), dimension(:,:), allocatable :: xout ! State output 
    end type IntegratorOut
    
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

    ! Main integrator procedure
    subroutine rk45(f_ptr, intin, intout)

        use kindmodule
        use mathmodule
        use derivkepmodule

        ! DECLARATION
        implicit none

        ! In/Out
        type(IntegratorIn),intent(in) :: intin 
        type(IntegratorOut),intent(out) :: intout

        abstract interface
        function func(t, x)
        use kindmodule
        complex(wp),intent(in) :: t
        complex(wp),intent(in),allocatable :: x(:)
        complex(wp),allocatable :: func(:)
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
    complex(wp),dimension(:),allocatable :: xtemp
    complex(wp),dimension(:),allocatable :: xdiff

    allocate(f(6,intin%n))
    allocate(x0(intin%n))
    allocate(xhatadd(intin%n))
    allocate(xadd(intin%n))
    allocate(xhat(intin%n))
    allocate(x(intin%n))
    allocate(xstore(intin%stepmax,intin%n))
    allocate(tstore(intin%stepmax))
    allocate(estore(intin%stepmax))
    allocate(hstore(intin%stepmax))
    allocate(xtemp(intin%n))
    allocate(xdiff(intin%n))

    ! EXECUTION

    ! Initialization
    x = intin%x0
    h = intin%h0
    step = 1

    do
        ! Check if we have reached final time or maximum steps
        if ((real(t) >= real(intin%tf)) .OR. (step > intin%stepmax)) then
            exit
        end if

        ! Evaluate derivatives
        f(1,:) = h*f_ptr(t, x)

        xtemp = x + beta(1,1)*f(1,:)
        f(2,:) = h*f_ptr(t + alpha(1)*h, xtemp)

        xtemp = x + beta(2,1)*f(1,:) + beta(2,2)*f(2,:)
        f(3,:) = h*f_ptr(t + alpha(2)*h, xtemp)
        
        xtemp = x + beta(3,1)*f(1,:) + beta(3,2)*f(2,:) + beta(3,3)*f(3,:)
        f(4,:) = h*f_ptr(t + alpha(3)*h, xtemp)
        
        xtemp = x + beta(4,1)*f(1,:) + beta(4,2)*f(2,:) + beta(4,3)*f(3,:) &
            + beta(4,4)*f(4,:) 
        f(5,:) = h*f_ptr(t + alpha(4)*h, xtemp)

        xtemp = x + beta(5,1)*f(1,:) + beta(5,2)*f(2,:) + beta(5,3)*f(3,:) &
            + beta(5,4)*f(4,:) + beta(5,5)*f(5,:)
        f(6,:) = h*f_ptr(t + alpha(5)*h, xtemp)

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
        xdiff = x - xhat
        error = astnorm(xdiff)

        ! Save current step
        xstore(step,:) = x
        tstore(step) = t
        estore(step) = error
        hstore(step) = h

        ! Step forward 
        if ((abs(real(error)) > intin%tol) .AND. (real(h) > intin%hmin)) then
            h = max(real(h)*intin%rparam, intin%hmin)
        else if (abs(real(error)) < (intin%tol)*intin%dbparam) then
            h = h*intin%iparam
        end if

        t = t + h
        step = step + 1

    end do 

    ! Populate output
    allocate(intout%xout(step-1,intin%n))
    allocate(intout%tout(step-1))
    allocate(intout%eout(step-1))
    allocate(intout%hout(step-1))

    do i=1,step-1
        intout%xout(i,:) = xstore(i,:)
        intout%tout(i) = tstore(i) 
        intout%eout(i) = estore(i)
        intout%hout(i) = hstore(i) 
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
    deallocate(xtemp)
    deallocate(xdiff)

end subroutine rk45

end module rk45module
