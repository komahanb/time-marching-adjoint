[![Build Status](https://travis-ci.org/komahanb/time-marching-adjoint.svg?branch=master)](https://travis-ci.org/komahanb/time-marching-adjoint)

# Time Marching and Unsteady Discrete Adjoint Sensitivities

Solves ordinary differential equations of the form: $R(\ldots, \ddot{q}, \dot{q}, q) = 0$

## How to build and run test case?

1. Within Makefile.in make sure F90 variable is set to use the fortran compiler that you have. gfortran > 6 is required.
2. cd src
3.  make
4 ./test_marching

## Implicit Time Marching Schemes:

 o Backward difference formulas
 o Adams Bashforth--Moulton
 o Diagonally Implicit Runge Kutta

## Unsteady Discrete Adjoint [unimplemented]

[Citation] (https://arc.aiaa.org/doi/abs/10.2514/6.2017-1671)