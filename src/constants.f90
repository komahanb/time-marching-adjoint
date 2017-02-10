#include "scalar.fpp"

!=====================================================================!
! Module defining all constants used
!=====================================================================!

module constants

  use iso_fortran_env , only : WP => real64

  implicit none

  ! Versioning numbers
  integer, parameter :: VERSION_MAJOR   = 0
  integer, parameter :: VERSION_MINOR   = 0
  integer, parameter :: VERSION_RELEASE = 0

  ! Physical constants
  type(scalar), parameter :: PI = 22.0_WP/7.0_WP !3.1415926535898_WP
  
  ! Handy double numbers
  type(scalar), parameter :: & 
       INFINITY     = huge(0.0_WP),       & ! biggest floating number
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
