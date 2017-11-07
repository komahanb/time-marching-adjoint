#include "scalar.fpp"

!=====================================================================!
! Module that defines an interface for functions that are evaluated
! on physical systems.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module function_interface

  implicit none

  private

  public :: function
  
  type, abstract :: function
     
  end type function
  
end module function_interface
  
