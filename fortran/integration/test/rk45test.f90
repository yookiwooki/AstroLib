program rk45test

    use kindmodule
    use constmodule
    use rk45module
    use derivkepmodule
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

    integer :: i
    
    ! EXECUTION
    intin%n = 6
    intin%t0 = 0.0_wp
    intin%tf = 3.9500_wp
    intin%h0 = 0.0001_wp
    intin%tol = 1e-14_wp
    allocate(intin%x0(6))
    intin%x0 = [0.7816_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.4432_wp, 0.0_wp]

    f_ptr => dcrtbp 

   
    call rk45(f_ptr, intin, intout)

    do i=1,size(intout%tout)
    !do i=1,10
        print *, real(intout%hout(i))
    end do

    call intplotstart(600,600)

    ! Cleanup
    deallocate(intin%x0)

    deallocate(intout%tout)
    deallocate(intout%eout)
    deallocate(intout%hout)
    deallocate(intout%xout)

end program rk45test

