#include "scalar.fpp"

!=====================================================================!
! Module defining all constants used
!=====================================================================!

module constants

  use iso_fortran_env

  implicit none

  ! Floating precision
  integer, parameter :: SP = REAL32  
  integer, parameter :: DP = REAL64
  integer, parameter :: QP = REAL128

  ! Set the working precision
  integer, parameter :: WP = DP

  ! Versioning numbers
  integer, parameter :: VERSION_MAJOR   = 0
  integer, parameter :: VERSION_MINOR   = 0
  integer, parameter :: VERSION_RELEASE = 0

  ! Physical constants
  type(scalar), parameter :: PI = 4.0_WP*atan(1.0_WP)
  !3.141592653589793238462643383279502884197169399375105820974944592307816406286
  !3.1415926535897931

  ! Handy double numbers
  type(scalar), parameter :: & 
       LARGE        = huge(0.0_WP),       & ! biggest floating number
       TINY         = epsilon(0.0_WP),    & ! smallest floating number
       ZERO         = 0.0_WP,             & ! zero
       ONE          = 1.0_WP,             & ! one
       TWO          = 2.0_WP,             & ! two
       HALF         = 0.5_WP                ! half

  ! Maximum number of words in a single line, length of line, and
  ! length of single word
  integer, parameter :: MAX_WORDS    = 500
  integer, parameter :: MAX_LINE_LEN = 250
  integer, parameter :: MAX_WORD_LEN = 150

  ! Unit numbers
  integer, parameter :: LOG_UNIT = 6 ! unit # for writing log file

end module constants
