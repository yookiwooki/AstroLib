program rk45test

    use kindmodule
    use constmodule
    use rk45module
    use derivkepmodule

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

    integer :: i
    
    ! EXECUTION
    intin%n = 6
    intin%t0 = 0.0
    intin%tf = 2.0*pi
    intin%h0 = 0.001
    allocate(intin%x0(6))
    intin%x0 = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]

    f_ptr => dkep

    call rk45(f_ptr, intin, intout)


    !do i=1,size(intout%tout)
    do i=1,10
        print *, real(intout%eout(i))
    end do


    ! Cleanup
    deallocate(intin%x0)

    deallocate(intout%tout)
    deallocate(intout%eout)
    deallocate(intout%hout)
    deallocate(intout%xout)

end program rk45test

