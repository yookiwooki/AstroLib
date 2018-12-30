!derivmodule.f90 - equations of motion for integrator
module derivmodule

    use kindmodule

    implicit none

    complex(wp) :: mu=1.0              ! Gravitational parameter
    complex(wp) :: mustar=0.0121505856 ! Mass ratio for three body problem 

    public :: dkep, dcrtbp, setmu, setmustar, getmu, getmustar
    private :: mu, mustar

contains

    ! Equations of motion for Keplerian dynamics
    subroutine dkep(t, x, xdot)

        use kindmodule
        use mathmodule

        ! DECLARATION
        implicit none

        ! Input/output
        complex(wp),intent(in) :: t
        complex(wp),intent(in),dimension(:) :: x
        complex(wp),intent(out),dimension(:) :: xdot
        ! Local
        complex(wp) :: rmag3

        ! EXECUTION
        xdot(1:3) = x(4:6)

        rmag3 = astnorm(x(1:3))**3
        xdot(4) = -mu*x(1)/rmag3
        xdot(5) = -mu*x(2)/rmag3
        xdot(6) = -mu*x(3)/rmag3

    end subroutine dkep


    ! Equations of motion for CRTBP 
    subroutine dcrtbp(t, x, xdot)

        use kindmodule
        use mathmodule

        ! DECLARATION
        implicit none

        ! Input/output
        complex(wp),intent(in) :: t
        complex(wp),intent(in),dimension(:) :: x
        complex(wp),intent(out),dimension(:) :: xdot
        ! Local
        complex(wp),dimension(3) :: r1, r2
        complex(wp) :: r1mag3, r2mag3

        ! EXECUTION
        xdot(1:3) = x(4:6)

        r1(1) = x(1) + mustar
        r1(2) = x(2)
        r1(3) = x(3)

        r2(1) = x(1) - 1 + mustar
        r2(2) = x(2)
        r2(3) = x(3)

        r1mag3 = astnorm(r1)**3
        r2mag3 = astnorm(r2)**3

        xdot(4) = 2*x(5) + x(1) - (1-mustar)*(x(1)+mustar)/r1mag3 - mustar*(x(1)-1+mustar)/r2mag3 
        xdot(5) = -2*x(4) + x(2) - (1-mustar)*x(2)/r1mag3 - mustar*x(2)/r2mag3 
        xdot(6) = -(1-mustar)*x(3)/r1mag3 - mustar*x(3)/r2mag3 

    end subroutine dcrtbp


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
