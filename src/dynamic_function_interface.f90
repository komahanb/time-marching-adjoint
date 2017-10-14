#include "scalar.fpp"

!=====================================================================!
! Module that defines an interface for time-dependent functions that
! are evaluated on physical systems.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module dynamic_function_interface

  use function_interface, only : function
  
  implicit none

  private

  public :: dynamic_function
  
  type, abstract, extends(function) :: dynamic_function
     
  end type dynamic_function
  
end module dynamic_function_interface
  
