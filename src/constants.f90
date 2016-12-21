#include "scalar.fpp"

!=====================================================================!
! Module defining all constants used
!=====================================================================!

module constants

  implicit none

  ! Versioning numbers
  integer, parameter :: VERSION_MAJOR = 0
  integer, parameter :: VERSION_MINOR = 0
  integer, parameter :: VERSION_RELEASE = 0

  ! Physical constants
  type(scalar), parameter :: PI = 3.1415926535898d0
  
  ! Handy double numbers
  type(scalar), parameter :: & 
       INFINITY     = huge(0.0d0),       & ! positive infinity
       ZERO         = 0.0d0,             & ! double zero
       ONE          = 1.0d0,             & ! double one
       TWO          = 2.0d0,             & ! double two
       HALF         = 0.5d0                ! double half

  ! Maximum number of words in a single line, length of line, and
  ! length of single word
  integer, parameter :: MAX_WORDS    = 500
  integer, parameter :: MAX_LINE_LEN = 250
  integer, parameter :: MAX_WORD_LEN = 150

  ! Unit numbers
  integer, parameter :: LOG_UNIT = 6 ! unit # for writing log file

end module constants
