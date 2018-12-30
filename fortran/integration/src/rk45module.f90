! ****************************************************************************************
! Author: Sean McArdle
! rk45module.f90 - Runge-Kutta-Fehlberg 4(5) integrator 

module rk45module

    use kindmodule

    implicit none

    private

    public :: IntegratorIn
    public :: IntegratorOut
    public :: rk45

    ! Input for integrator 
    type IntegratorIn
        integer :: n                               ! State size
        complex(wp) :: t0                          ! Initial time
        complex(wp) :: tf                          ! Final time
        complex(wp) :: h0=1.0E-8_wp                ! Initial step size
        integer :: stepmax = 100000                ! Max number of steps
        real(wp) :: hmin = 1.0E-10_wp               ! Min step size
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
! ****************************************************************************************

! ****************************************************************************************
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
            subroutine func(t,x,xdot)
                use kindmodule
                complex(wp),intent(in) :: t
                complex(wp),intent(in),dimension(:) :: x
                complex(wp),intent(out),dimension(:) :: xdot 
            end subroutine func
        end interface 
        procedure (func), pointer :: f_ptr => null ()  ! Derivative function 

        ! Local variables
        complex(wp) :: t                                 ! Time
        complex(wp),dimension(:),allocatable :: xdot     ! Derivatives
        complex(wp),dimension(:),allocatable :: f1,f2,f3 ! Derivatives
        complex(wp),dimension(:),allocatable :: f4,f5,f6 ! Derivatives
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

        allocate(xdot(intin%n))
        allocate(f1(intin%n))
        allocate(f2(intin%n))
        allocate(f3(intin%n))
        allocate(f4(intin%n))
        allocate(f5(intin%n))
        allocate(f6(intin%n))
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
            call f_ptr(t, x, xdot)
            f1 = h*xdot

            xtemp = x + beta(1,1)*f1
            call f_ptr(t + alpha(1)*h, xtemp, xdot)
            f2 = h*xdot

            xtemp = x + beta(2,1)*f1 + beta(2,2)*f2
            call f_ptr(t + alpha(2)*h, xtemp, xdot)
            f3 = h*xdot

            xtemp = x + beta(3,1)*f1 + beta(3,2)*f2 + beta(3,3)*f3
            call f_ptr(t + alpha(3)*h, xtemp, xdot)
            f4 = h*xdot

            xtemp = x + beta(4,1)*f1 + beta(4,2)*f2 + beta(4,3)*f3 &
                + beta(4,4)*f4 
            call f_ptr(t + alpha(4)*h, xtemp, xdot)
            f5 = h*xdot

            xtemp = x + beta(5,1)*f1 + beta(5,2)*f2 + beta(5,3)*f3 &
                + beta(5,4)*f4 + beta(5,5)*f5
            call f_ptr(t + alpha(5)*h, xtemp, xdot)
            f6 = h*xdot

            ! Evaluate integrated states
            xhat = x + chat(1)*f1 + chat(2)*f2 + chat(3)*f3 &
                + chat(4)*f4 + chat(5)*f5 + chat(6)*f6

            x = x + c(1)*f1 + c(2)*f2 + c(3)*f3 &
                + c(4)*f4 + c(5)*f5
            
            ! Check difference between high/low order 
            xdiff = x - xhat
            error = astnorm(xdiff)/astnorm(x)

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
        deallocate(xdot)
        deallocate(f1)
        deallocate(f2)
        deallocate(f3)
        deallocate(f4)
        deallocate(f5)
        deallocate(f6)
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
! ****************************************************************************************
    end module rk45module
! ****************************************************************************************
