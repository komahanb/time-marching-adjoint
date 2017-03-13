#include "scalar.fpp"

!=====================================================================!
! Module that contains common procedures for any physical system
! subject to governing equations
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module physics_interface

  implicit none
  
  private

  public :: physics
 
  !-------------------------------------------------------------------!
  ! Type that models any physical system
  !-------------------------------------------------------------------!
  
  type, abstract :: physics

     type(integer) :: num_state_vars 

   contains  

     ! Provided procedures
     procedure :: get_num_state_vars
     procedure :: set_num_state_vars
     
     ! Deferred procedures
     procedure(add_residual_interface), deferred :: add_residual
     procedure(add_jacobian_interface), deferred :: add_jacobian

  end type physics

  ! Interfaces to deferred procedures
  abstract interface

     !----------------------------------------------------------------!
     ! Interface for residual assembly
     !----------------------------------------------------------------!

     pure subroutine add_residual_interface(this, residual)
       
       import :: physics
       
       class(physics), intent(inout) :: this
       type(scalar)  , intent(inout) :: residual(:)
       
     end subroutine add_residual_interface

     !----------------------------------------------------------------!
     ! Interface for jacobian assembly
     !----------------------------------------------------------------!
     
     pure subroutine add_jacobian_interface(this, jacobian, coeff)

       import :: physics

       class(physics), intent(inout) :: this
       type(scalar)  , intent(inout) :: jacobian(:,:)
       type(scalar)  , intent(in)    :: coeff(:)

     end subroutine add_jacobian_interface

  end interface

contains
 
  !===================================================================!
  ! Returns the number of state variables in the physical system
  !===================================================================!
  
  pure type(integer) function get_num_state_vars(this)

    class(physics), intent(in) :: this

    get_num_state_vars = this % num_state_vars

  end function get_num_state_vars

  !===================================================================!
  ! Sets the number of state variables in the physical system
  !===================================================================!
  
  pure subroutine set_num_state_vars(this, num_state_vars)

    class(physics), intent(inout) :: this
    type(integer)  , intent(in)    :: num_state_vars

    this % num_state_vars  = num_state_vars

  end subroutine set_num_state_vars

end module physics_interface
