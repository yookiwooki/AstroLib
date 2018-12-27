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

        ! EXECUTION
        astnorm = sqrt(dot_product(x,x))

    end function astnorm

end module mathmodule
