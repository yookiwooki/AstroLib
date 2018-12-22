! kindmodule.f90 
! Defines working precision

module kindmodule

    use, intrinsic :: iso_fortran_env

    implicit none

    integer, parameter, public :: wp=real64
    integer, parameter, public :: sp=kind(1e0)

end module kindmodule
