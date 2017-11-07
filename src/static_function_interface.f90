#include "scalar.fpp"

!=====================================================================!
! Module that defines an interface for time-independent functions that
! are evaluated on physical systems.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module static_function_interface

  use function_interface, only : function
  
  implicit none

  private

  public :: static_function
  
  type, abstract, extends(function) :: static_function
     
  end type static_function
  
end module static_function_interface
  
