#include "scalar.fpp"

!=====================================================================!
! Module that contains common procedures for any time-independent
! physical system subject to governing equations of n-th order
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module time_independent_physics_interface

  use constants, only : WP
  use physics_interface, only : physics

  implicit none
  
  private
  
  public :: statics
  
  type, abstract, extends(physics) :: statics
    
   contains

  end type statics
  
contains
  
end module time_independent_physics_interface
