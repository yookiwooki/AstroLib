! intplotmodule.f90 - plotting utility for integrator
module intplotmodule

    use opengl_gl
    
    implicit none

    integer :: cpoint=0 
    real(glfloat),dimension(:) :: X(10000), Y(10000)
    real(glfloat) :: intscale=0.9

contains

    subroutine addpoint(xpoint,ypoint)
        use opengl_gl

        ! DECLARATION
        implicit none

        real(glfloat), intent(in) :: xpoint, ypoint

        ! EXECUTION
        cpoint = cpoint + 1
        X(cpoint) = xpoint
        Y(cpoint) = ypoint

    end subroutine addpoint


    subroutine drawcircle(xcenter, radius)
        use opengl_gl
        use opengl_glut

        ! DECLARATION
        implicit none 

        integer(glint) :: i, n
        real(glfloat), intent(in) :: xcenter, radius
        real(glfloat) :: xcirc, ycirc, xcoord, ycoord
        real(glfloat),parameter :: twicepi = 3.14159*2.0

        !call glclear(GL_COLOR_BUFFER_BIT + GL_DEPTH_BUFFER_BIT)
        call glcolor3f(1.0_glfloat, 0.0_glfloat, 0.0_glfloat)

        xcirc = xcenter 
        ycirc = 0.0
        n = 20

        ! Draw planet
        call glbegin(GL_TRIANGLE_FAN)
        call glvertex2f(xcirc*intscale, ycirc*intscale)
        do i = 1, n+1
            xcoord = xcirc + radius*cos((real(i)*twicepi)/real(n))
            ycoord = ycirc + radius*sin((real(i)*twicepi)/real(n))
            call glvertex2f(xcoord*intscale, ycoord*intscale)
        end do
        call glend
        call glflush

    end subroutine drawcircle 


    subroutine drawtraj()
        use opengl_gl
        use opengl_glut

        ! DECLARATION
        implicit none 

        integer :: i

        ! Draw trajectory
        call glcolor3f(0.0_glfloat, 0.0_glfloat, 1.0_glfloat)
        call glbegin(GL_LINE_STRIP)
        do i=1,cpoint 
            call glvertex2f(X(i)*intscale,Y(i)*intscale)
        end do
        call glend
        call glflush
    end subroutine drawtraj 


    subroutine mainloop()
        use opengl_gl
        use opengl_glut

        ! DECLARATION
        implicit none

        ! EXECUTION
        call glloadidentity
        call gltranslatef(-1.0_glfloat*intscale, 0.0_glfloat, 0.0_glfloat)
        call drawcircle(0.0, 0.02)
        call drawcircle(1.0, 0.01)
        call drawtraj
        call glutswapbuffers

    end subroutine mainloop


    subroutine glsetup(width, height)

        use opengl_gl
        use opengl_glu
        use opengl_glut
        
        ! DECLARATION
        implicit none

        integer(glint), intent(in) :: width, height
        real(gldouble) :: aspect

        ! EXECUTION
        aspect = real(width)/real(height)
        call glviewport(-1, -1, width, height)
        call glmatrixmode(GL_PROJECTION)
        call glenable(GL_DEPTH_TEST)
        call gluperspective(45.0_gldouble, aspect, 0.1_gldouble, 100.0_gldouble)

    end subroutine glsetup


    subroutine intplotstart(width, height)
        
        use opengl_gl
        use opengl_glut

        ! DECLARATION
        implicit none

        integer(glint), intent(in) :: width, height
        integer :: i

        ! EXECUTION
        call glutinit()
        call glutinitwindowsize(width, height)
        call glutinitdisplaymode(GLUT_RGB + GLUT_DOUBLE)
        i = glutcreatewindow("Trajectory")
        call glutdisplayfunc(mainloop)
        call glsetup(width, height)
        call glutmainloop()

    end subroutine intplotstart


end module intplotmodule

