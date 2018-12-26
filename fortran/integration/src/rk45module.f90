! rk45module.f90 - Runge-Kutta-Fehlberg 4(5) integrator 
module rk45module

    use kindmodule
    use intplotmodule

    implicit none

    ! Input for integrator 
    type IntegratorIn
        integer :: n                               ! State size
        complex(wp) :: t0                          ! Initial time
        complex(wp) :: tf                          ! Final time
        complex(wp) :: h0=1.0E-3_wp                ! Initial step size
        integer :: stepmax = 100000                ! Max number of steps
        real(wp) :: hmin = 1.0E-6_wp               ! Min step size
        real(wp) :: tol = 1.0E-12_wp               ! Tolerance
        complex(wp),dimension(:),allocatable :: x0 ! Initial state 
        real(wp) :: rparam = 0.5_wp                ! Step reduce factor
        real(wp) :: iparam = 1.01_wp               ! Step increase factor
        real(wp) :: dbparam = 0.8_wp               ! Step increase deadband
    end type IntegratorIn

    ! Output for integrator
    type IntegratorOut
        complex(wp),dimension(:),allocatable :: tout   ! Time
        complex(wp),dimension(:),allocatable :: eout   ! Error
        complex(wp),dimension(:),allocatable :: hout   ! Step size
        complex(wp),dimension(:,:),allocatable :: xout ! State
    end type IntegratorOut
    
    ! Runge-Kutta-Fehlberg coefficients (see reference)
    complex(wp), parameter :: alpha(5) = &
        [1.0_wp/4.0_wp, 3.0_wp/8.0_wp, 12.0_wp/13.0_wp, 1.0_wp, 1.0_wp/2.0_wp]

    complex(wp), parameter :: beta(5,5) = reshape([ &
        1.0_wp/4.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
        3.0_wp/32.0_wp, 9.0_wp/32.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
        1932.0_wp/2197.0_wp, -7200.0_wp/2197.0_wp, 7296.0_wp/2197.0_wp, 0.0_wp, 0.0_wp, &
        439.0_wp/216.0_wp, -8.0_wp, 3680.0_wp/513.0_wp, -845.0_wp/4104.0_wp, 0.0_wp, &
        -8.0_wp/27.0, 2.0_wp, -3544.0_wp/2565.0_wp, 1859.0_wp/4104.0_wp, -11.0_wp/40.0_wp], &
        shape(beta),order=[2,1])

    complex(wp), parameter :: c(5) = &
        [25.0_wp/216.0_wp, 0.0_wp, 1408.0_wp/2565.0_wp, 2197.0_wp/4104.0_wp, -1.0_wp/5.0_wp]

    complex(wp), parameter :: chat(6) = &
        [16.0_wp/135.0_wp, 0.0_wp, 6656.0_wp/12825.0_wp, 28561.0_wp/56430.0_wp, &
        -9.0_wp/50.0_wp, 2.0_wp/55.0_wp]

contains 

    ! RKF 4(5)th order variable step size numerical integration 
    ! 
    ! Ref:
    ! Fehlberg, E. "Low order classical runge-kutta formulas with stepsize  
    ! control and their application to some heat transfer problems." (1969).
    ! 
    ! NOTE: The reference above has a tricky typo related to the
    ! step size in Equation 2 that has been corrected in this implementation.

    subroutine rk45(f_ptr, intin, intout)

        use kindmodule
        use mathmodule
        use derivmodule

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

        ! Local variables
        complex(wp) :: t                                 ! Time
        complex(wp),dimension(:,:),allocatable :: f      ! Derivatives
        complex(wp) :: h                                 ! Step size
        complex(wp),dimension(:),allocatable :: xhat     ! High order result
        complex(wp),dimension(:),allocatable :: x        ! Low order result
        complex(wp),dimension(:),allocatable :: xold     ! Previous iter result
        complex(wp),dimension(:,:),allocatable :: xstore ! State storage
        complex(wp),dimension(:),allocatable :: tstore   ! Time storage
        complex(wp),dimension(:),allocatable :: estore   ! Error storage
        complex(wp),dimension(:),allocatable :: hstore   ! Step size storage
        integer :: step                                  ! Step counter
        complex(wp) :: error                             ! High/low order error
        integer :: i                                     ! Loop index
        complex(wp),dimension(:),allocatable :: xtemp    ! State input for deriv
        complex(wp),dimension(:),allocatable :: xdiff    ! High/low order diff

        allocate(f(6,intin%n))
        allocate(xhat(intin%n))
        allocate(x(intin%n))
        allocate(xold(intin%n))
        allocate(xstore(intin%stepmax,intin%n)) ! Store output in temporary
        allocate(tstore(intin%stepmax))         ! arrays with the maximum
        allocate(estore(intin%stepmax))         ! possible elements.
        allocate(hstore(intin%stepmax))
        allocate(xtemp(intin%n))
        allocate(xdiff(intin%n))

        ! EXECUTION

        ! Initialization
        x = intin%x0
        t = intin%t0
        h = intin%h0
        step = 1

        ! Plot first step
        call resetpoints
        call addpoint(real(x(1), sp), real(x(2), sp))

        ! Save first step
        xstore(step,:) = x
        tstore(step) = t
        estore(step) = 0 
        hstore(step) = 0
        step = step + 1

        do
            t = t + h
            xold = x

            ! Check if we have reached final time or maximum steps
            if (real(t) > real(intin%tf)) then
                exit
            elseif (step > intin%stepmax) then
                print *, "WARNING: stepmax reached (rk45module)"
                exit
            elseif (real(t+h) > real(intin%tf)) then
                t = intin%tf 
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
            xhat = x + chat(1)*f(1,:) + chat(2)*f(2,:) + chat(3)*f(3,:) &
                + chat(4)*f(4,:) + chat(5)*f(5,:) + chat(6)*f(6,:)

            x = x + c(1)*f(1,:) + c(2)*f(2,:) + c(3)*f(3,:) &
                + c(4)*f(4,:) + c(5)*f(5,:)
            
            ! Check difference between high/low order 
            xdiff = x - xhat
            error = astnorm(xdiff)

            ! Plot current step
            call addpoint(real(x(1), sp), real(x(2), sp))

            ! Save current step
            xstore(step,:) = x
            tstore(step) = t
            estore(step) = error
            hstore(step) = h

            ! Decrease step size if error above tolerance and above min step size
            if ((abs(real(error)) > intin%tol) .AND. (real(h) > intin%hmin)) then
                t = t - h
                x = xold
                if (real(h)*intin%rparam <= intin%hmin) then
                    print *, "WARNING: minimum step size reached (rk45module)"
                end if
                h = max(real(h)*intin%rparam, intin%hmin)
            ! Increase step size if error is far enough above tolerance
            else if (abs(real(error)) < (intin%tol)*intin%dbparam) then
                h = h*intin%iparam
                step = step + 1
            ! Otherwise continue without changing step size 
            else
                step = step + 1
            end if

        end do 

        ! Populate output
        allocate(intout%xout(step-1,intin%n))
        allocate(intout%tout(step-1))
        allocate(intout%eout(step-1))
        allocate(intout%hout(step-1))

        intout%xout(:,:) = xstore(1:step-1,:)
        intout%tout(:) = tstore(1:step-1) 
        intout%eout(:) = estore(1:step-1)
        intout%hout(:) = hstore(1:step-1) 

        ! Cleanup
        deallocate(f)
        deallocate(xhat)
        deallocate(x)
        deallocate(xold)
        deallocate(xstore)
        deallocate(tstore)
        deallocate(estore)
        deallocate(hstore)
        deallocate(xtemp)
        deallocate(xdiff)

    end subroutine rk45

    end module rk45module
