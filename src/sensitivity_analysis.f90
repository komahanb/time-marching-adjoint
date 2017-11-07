#include "scalar.fpp"

!=====================================================================!
! Module that defines an interface for SENSITIVITY ANALYSIS of
! physical systems. This interface is a 'decorator' on top of analysis
! interface.
!
! Main steps in sensitivty analysis are:
!
! (a) solving the physical system for state variables
! (b) solving the adjoint/direct sensitivity variables
! 
! (i)  evaluating the functions of interest
! (ii) evaluating the derivatives of functions of interest
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module sensitivity_analysis_interface

  use analysis_interface, only : ianalysis => analysis
  use function_interface, only : ifunction => function

  implicit none
  
  private
  
  public :: sensitivity_analysis
  
  !-------------------------------------------------------------------!
  ! Type definition for analysis
  !-------------------------------------------------------------------!
  
  type, extends(analysis), abstract :: sensitivity_analysis

     ! Sentivity analysis needs the behavior of 'analysis' class
     class(ianalysis), allocatable :: analysis ! used to evaluate the states
     class(ifunction), allocatable :: function ! used to evaluate functions of interest

   contains  

     ! Implement methods to evalaute the sensitivity variables
     
  end type sensitivity_analysis

  ! Interfaces to deferred procedures
  abstract interface
  end interface

contains
  
end module sensitivity_analysis_interface
