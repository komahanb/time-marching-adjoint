#include "scalar.fpp"

!=====================================================================!
! Module that contains common procedures for any time-dependent
! physical system subject to governing equations of n-th order
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module dynamic_analysis_interface

  use analysis_interface, only : analysis
  use dynamic_physics_interface, only : dynamics

  implicit none
  
  private
  
  public :: dynamic_analysis
  
  type, abstract, extends(analysis) :: dynamic_analysis
     
   contains
     
  end type dynamic_analysis

  ! Interfaces to deferred procedures
  interface

  end interface

contains
  
end module dynamic_analysis_interface
