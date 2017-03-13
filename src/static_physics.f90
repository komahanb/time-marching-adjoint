#include "scalar.fpp"

!=====================================================================!
! Module that contains common procedures for any time-independent
! physical system subject to governing equations of n-th order
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module static_physics_interface

  use physics_interface, only : physics

  implicit none
  
  private
  
  public :: statics
  
  type, abstract, extends(physics) :: statics
    
   contains

  end type statics
  
contains
  
end module static_physics_interface
