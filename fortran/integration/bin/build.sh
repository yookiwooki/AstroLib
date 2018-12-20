#!/usr/bin/env bash

gfortran -g -O0 ../../kindmodule.f90 ../../constmodule.f90 ../../mathmodule.f90 \
    ../src/derivkepmodule.f90 ../src/rk45module.f90 \
    ../test/rk45test.f90 -o rk45test

