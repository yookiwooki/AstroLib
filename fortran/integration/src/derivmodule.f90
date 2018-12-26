!derivmodule.f90 - equations of motion for integrator
module derivmodule
    use kindmodule

    implicit none

    complex(wp) :: mu=1.0              ! Gravitational parameter
    complex(wp) :: mustar=0.0121505856 ! Mass ratio for three body problem 

    PUBLIC :: dkep, dcrtbp, setmu, setmustar, getmu, getmustar
    PRIVATE :: mu, mustar

contains

    ! Equations of motion for Keplerian dynamics
    function dkep(t, x)

        use kindmodule
        use mathmodule

        ! DECLARATION
        implicit none

        complex(wp),intent(in) :: t
        complex(wp),intent(in),allocatable :: x(:)
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


    ! Equations of motion for CRTBP 
    function dcrtbp(t, x)

        use kindmodule
        use mathmodule

        ! DECLARATION
        implicit none

        complex(wp),intent(in) :: t
        complex(wp),intent(in),allocatable :: x(:)
        complex(wp),dimension(3) :: accel
        complex(wp),allocatable :: dcrtbp(:)
        complex(wp),allocatable :: r1(:), r2(:)
        complex(wp) :: r1mag3, r2mag3

        ! EXECUTION
        if (.NOT. allocated(r1)) then 
            allocate(r1(3))
        end if

        if (.NOT. allocated(r2)) then 
            allocate(r2(3))
        end if

        if (.NOT. allocated(dcrtbp)) then 
            allocate(dcrtbp(6))
        end if

        dcrtbp(1:3) = x(4:6)

        r1(1) = x(1) + mustar
        r1(2) = x(2)
        r1(3) = x(3)

        r2(1) = x(1) - 1 + mustar
        r2(2) = x(2)
        r2(3) = x(3)

        r1mag3 = astnorm(r1)**3
        r2mag3 = astnorm(r2)**3

        accel(1) = 2*x(5) + x(1) - (1-mustar)*(x(1)+mustar)/r1mag3 - mustar*(x(1)-1+mustar)/r2mag3 
        accel(2) = -2*x(4) + x(2) - (1-mustar)*x(2)/r1mag3 - mustar*x(2)/r2mag3 
        accel(3) = -(1-mustar)*x(3)/r1mag3 - mustar*x(3)/r2mag3 
        dcrtbp(4:6) = accel

    end function dcrtbp


    ! Sets mass ratio module variable
    subroutine setmustar(mustarin)
        use kindmodule
        implicit none
        complex(wp),intent(in) :: mustarin
        mustar = mustarin
    end subroutine setmustar


    ! Reports mass ratio module variable
    function getmustar()
        use kindmodule
        implicit none
        complex(wp) :: getmustar
        getmustar = mustar
    end function getmustar


    ! Sets gravitational parameter module variable
    subroutine setmu(muin)
        use kindmodule
        implicit none
        complex(wp),intent(in) :: muin
        mu = muin
    end subroutine setmu


    ! Reports gravitational parameter module variable
    function getmu()
        use kindmodule
        implicit none
        complex(wp) :: getmu
        getmu = mu
    end function getmu


end module derivmodule
