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
    intin%t0 = 0.0_wp
    intin%tf = 18.0_wp
    intin%h0 = 0.001_wp
    allocate(intin%x0(6))
    intin%x0 = [1.0_wp, 0.0_wp, 0.0_wp, 0.1_wp, 1.0_wp, 0.1_wp]

    f_ptr => dkep

    call rk45(f_ptr, intin, intout)

    do i=1,size(intout%tout)
    !do i=1,10
        print *, real(intout%tout(i))
    end do

    ! Cleanup
    deallocate(intin%x0)

    deallocate(intout%tout)
    deallocate(intout%eout)
    deallocate(intout%hout)
    deallocate(intout%xout)

end program rk45test

