#include "scalar.fpp"

!=====================================================================!
! Module that contains common procedures for any time-dependent
! physical system subject to governing equations of n-th order
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module time_dependent_physics_interface

  use constants, only : WP
  use vector_interface, only : vector
  use physics_interface, only : physics

  implicit none
  
  private
  
  public :: dynamics
  
  type, abstract, extends(physics) :: dynamics

     type(integer) :: time_order ! order of the diff. eqn.
     
   contains

     procedure(initial_condition_interface), deferred :: get_initial_condition

  end type dynamics

  ! Interfaces to deferred procedures
  interface

     !----------------------------------------------------------------!
     ! Supplying the initial condition to march in time
     !----------------------------------------------------------------!
     
     pure subroutine initial_condition_interface(this, state_vectors)

       import :: dynamics, vector

       class(dynamics), intent(inout) :: this
       class(vector), intent(inout) :: state_vectors(:)

     end subroutine initial_condition_interface

  end interface

contains

end module time_dependent_physics_interface
