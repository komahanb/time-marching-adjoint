#include "scalar.fpp"

!=====================================================================!
! Static analysis of time-independent systems
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module static_analysis_interface

  use analysis_interface, only : analysis

  implicit none
  
  private
  
  public :: statics
  
  type, abstract, extends(analysis) :: statics
    
   contains

  end type statics
  
contains
  
end module static_analysis_interface
