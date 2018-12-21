module derivkepmodule

    implicit none

contains

    ! Equations of motion for Keplerian dynamics
    function dkep(t, x)

        use kindmodule
        use mathmodule

        ! DECLARATION
        implicit none

        complex(wp),intent(in) :: t
        complex(wp),intent(in),allocatable :: x(:)
        complex(wp),parameter :: mu=1.0
        complex(wp),dimension(3) :: accel
        complex(wp),allocatable :: dkep(:)
        complex(wp),allocatable :: pos(:)
        complex(wp) :: rmag3

        ! EXECUTION
        if (.NOT. allocated(pos)) then 
            allocate(pos(3))
        end if

        if (.NOT. allocated(dkep)) then 
            allocate(dkep(6))
        end if

        dkep(1:3) = x(4:6)

        pos(1:3) = x(1:3)
        rmag3 = astnorm(pos)**3
        accel(1) = -mu*x(1)/rmag3
        accel(2) = -mu*x(2)/rmag3
        accel(3) = -mu*x(3)/rmag3
        dkep(4:6) = accel

    end function dkep

end module derivkepmodule
