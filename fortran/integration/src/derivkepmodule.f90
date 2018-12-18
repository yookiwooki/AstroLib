module derivkepmodule

implicit none

contains

    ! Equations of motion for Keplerian dynamics
    function dkep(t, x)

        use kindmodule

        ! DECLARATION
        implicit none

        complex(wp),intent(in) :: t
        complex(wp),dimension(6),intent(in) :: x
        complex(wp),parameter :: mu=1.0
        complex(wp),dimension(3) :: accel
        complex(wp),dimension(6) :: dkep

        ! EXECUTION
        dkep(1:3) = x(4:6)
        accel(1) = -mu*x(1)/(norm2(real(x(1:3))**3))
        accel(2) = -mu*x(2)/(norm2(real(x(1:3))**3))
        accel(3) = -mu*x(3)/(norm2(real(x(1:3))**3))
        dkep(4:6) = accel

    end function dkep

end module derivkepmodule
