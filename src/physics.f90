#include "scalar.fpp"
!=====================================================================!
! Module that contains common procedures for any physical system
! subject to governing equations
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module physics_interface

  use constants, only : WP
  use vector_interface, only : vector
  use matrix_interface, only : matrix

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
     procedure(residual_assembly_interface), deferred :: assemble_residual
     procedure(jacobian_assembly_interface), deferred :: assemble_jacobian

  end type physics

  ! Interfaces to deferred procedures
  interface

     !----------------------------------------------------------------!
     ! Interface for residual assembly at each time step
     !----------------------------------------------------------------!

     pure subroutine residual_assembly_interface(this, residual, state_vectors)
       
       import :: physics, vector
       
       class(physics), intent(inout) :: this
       class(vector),  intent(inout) :: residual
       class(vector),  intent(in)    :: state_vectors(:)
       
     end subroutine residual_assembly_interface

     !----------------------------------------------------------------!
     ! Interface for jacobian assembly at each time step
     !----------------------------------------------------------------!
     
     pure subroutine jacobian_assembly_interface(this, jacobian, state_vectors, coeffs)

       import :: physics, vector, matrix

       class(physics) , intent(inout) :: this
       class(matrix)  , intent(inout) :: jacobian
       class(vector)  , intent(in)    :: state_vectors(:)
       type(scalar)   , intent(in)    :: coeffs(:)

     end subroutine jacobian_assembly_interface

  end interface

contains

end module physics_interface
