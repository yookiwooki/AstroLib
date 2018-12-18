program rk45test

    use kindmodule
    use rk45module
    use derivkepmodule

    ! DECLARATION
    implicit none

    abstract interface
    function func(t, x)
    use kindmodule
    complex(wp),intent(in) :: t
    complex(wp),dimension(6),intent(in) :: x
    complex(wp),dimension(6) :: func 
    end function func
    end interface 

    procedure (func), pointer :: f_ptr => null()
    type(IntegratorInOut) :: io

    integer :: i
    
    ! EXECUTION
    io%n = 6
    io%t0 = 0.0
    io%tf = 5.0
    allocate(io%x0(6))
    io%x0 = [1.0, 0.0, 0.0, 0.1, 1.0, 0.1]

    f_ptr => dkep

    call rk45(f_ptr, io)


    do i=1,size(io%tout)
        print *, real(io%xout(i,:))
    end do


    ! Cleanup
    deallocate(io%x0)
    deallocate(io%tout)
    deallocate(io%eout)
    deallocate(io%hout)
    deallocate(io%xout)

end program rk45test

