#include "scalar.fpp"

!=====================================================================!
! Module that contains common procedures for ANALYSIS of physical
! systems. The analysis is integration of the differential equations
! for time dependent systems and linear/nonlinear solve for
! time-independent systems.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module analysis_interface

  implicit none
  
  private

  public :: analysis
 
  !-------------------------------------------------------------------!
  ! Type definition for analysis
  !-------------------------------------------------------------------!
  
  type, abstract :: analysis

   contains  

  end type analysis

  ! Interfaces to deferred procedures
  abstract interface

  end interface

contains
  
end module analysis_interface
