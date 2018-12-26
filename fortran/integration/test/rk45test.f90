program rk45test

    use kindmodule
    use mathmodule
    use constmodule
    use rk45module
    use derivmodule
    use intplotmodule

    ! DECLARATION
    implicit none

    abstract interface
    function func(t, x)
    use kindmodule
    complex(wp),intent(in) :: t
    complex(wp),intent(in),allocatable :: x(:)
    complex(wp),allocatable :: func(:)
    end function func
    end interface 

    procedure (func), pointer :: f_ptr => null()
    type(IntegratorIn) :: intin 
    type(IntegratorOut) :: intout
    complex(wp),allocatable :: error(:)

    integer :: i
    
    ! EXECUTION
    intin%n = 6
    intin%t0 = 0.0_wp
    intin%tf = 2.0_wp*pi 
    intin%h0 = 1.0e-10_wp 
    intin%tol = 1e-15_wp
    allocate(intin%x0(6))
    !intin%x0 = [0.7816_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.4432_wp, 0.0_wp]
    intin%x0 = [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp]

    f_ptr => dkep
   
    call rk45(f_ptr, intin, intout)

    !do i=1,size(intout%tout)
    !!do i=1,10
    !    print *, real(intout%eout(i))
    !end do

    allocate(error(6))
    error = intin%x0 - intout%xout(size(intout%tout),:)
    print *, real(astnorm(error))
    print *, real(intin%x0(1:3))
    print *, real(intout%xout(size(intout%tout),1:3))

    call intplotstart(600,600)

    ! Cleanup
    deallocate(intin%x0)

    deallocate(intout%tout)
    deallocate(intout%eout)
    deallocate(intout%hout)
    deallocate(intout%xout)

end program rk45test

