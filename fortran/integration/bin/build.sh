#!/usr/bin/env bash

gfortran -g  ../../kindmodule.f90 ../../constmodule.f90 ../../mathmodule.f90 \
    ../plot/intplotmodule.f90 \
    ../src/derivkepmodule.f90 ../src/rk45module.f90 \
    -I$F03GL ${F03GL}OpenGL_gl.o ${F03GL}OpenGL_glu.o ${F03GL}OpenGL_freeglut.o \
    ../test/rk45test.f90 \
    -L/usr/lib/x86_64-linux-gnu -lglut -lGL -lGLU \
    -L/usr/lib/x86_64-linux-gnu -lXext -lX11 -lm \
    -o rk45test

