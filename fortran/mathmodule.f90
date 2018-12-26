! mathmodule.f90 - useful math functions for astrolib
module mathmodule

    implicit none

    contains

    ! Two-norm that works with complex numbers
    function astnorm(x)

        use kindmodule

        ! DECLARATION
        implicit none

        complex(wp),dimension(:),intent(in),allocatable :: x
        complex(wp) :: astnorm
        complex(wp) :: adder
        integer :: i

        ! EXECUTION
        adder = 0.0_wp
        do i=1,size(x)
            adder = adder + x(i)**2 
        end do
        astnorm = sqrt(adder)

    end function astnorm

end module mathmodule
