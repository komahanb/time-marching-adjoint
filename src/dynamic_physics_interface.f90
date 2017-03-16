#include "scalar.fpp"

!=====================================================================!
! Module that contains common procedures for any time-dependent
! physical system subject to governing equations of n-th order
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module dynamic_physics_interface

  use physics_interface, only : physics

  implicit none
  
  private
  
  public :: dynamics
  
  type, abstract, extends(physics) :: dynamics

     type(integer) :: time_deriv_order ! order of the differential equation
     
   contains
     
     ! Defined procedures
     procedure :: get_time_deriv_order
     procedure :: set_time_deriv_order

     ! Deferred procedure to subtypes
     procedure(initial_condition_interface), deferred :: get_initial_condition

  end type dynamics

  ! Interfaces to deferred procedures
  interface

     !----------------------------------------------------------------!
     ! Supplying the initial condition to march in time
     !----------------------------------------------------------------!
     
     pure subroutine initial_condition_interface(this, U)

       import :: dynamics

       class(dynamics), intent(in)  :: this
       type(scalar)   , intent(inout) :: U(:,:)

     end subroutine initial_condition_interface

  end interface

contains
  
  !===================================================================!
  ! Returns the highest order of time derivative in the physics
  !===================================================================!
  
  pure type(integer) function get_time_deriv_order(this)

    class(dynamics), intent(in) :: this

    get_time_deriv_order = this % time_deriv_order

  end function get_time_deriv_order

  !===================================================================!
  ! Sets the highest order of time derivative in the physics
  !===================================================================!
  
  pure subroutine set_time_deriv_order(this, order)

    class(dynamics), intent(inout) :: this
    type(integer), intent(in)      :: order

    this % time_deriv_order = order

  end subroutine set_time_deriv_order
  
end module dynamic_physics_interface
