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
     
     pure subroutine add_jacobian_interface(this, jacobian)

       import :: physics

       class(physics), intent(inout) :: this
       type(scalar)  , intent(inout) :: jacobian(:,:)

     end subroutine add_jacobian_interface

  end interface

contains

end module physics_interface
