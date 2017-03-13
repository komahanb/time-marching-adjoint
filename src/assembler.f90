!=====================================================================!
! Assemble the residual and jacobian of physical systems. This module
! is a wrapper on physical systems. The purpose of this module is to
! enclose the residual and jacobian assembly logic. The design is such
! that it will use a list of physical systems to assemble the global
! jacobian. The aim is to use for multidisciplinary solution
! approaches.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module assembler_interface

  use physics_interface , only : physics
  use matrix_interface  , only : matrix
  use vector_interface  , only : vector

  implicit none

  private

  public :: assembler

  type, abstract :: assembler

     class(physics), allocatable :: systems(:)

   contains

     ! Deferred procedures
     procedure(assemble_residual_interface), deferred :: assemble_residual
     procedure(assemble_jacobian_interface), deferred :: assemble_jacobian

  end type assembler

  ! Interfaces to deferred procedures
  abstract interface

     !----------------------------------------------------------------!
     ! Interface for residual assembly
     !----------------------------------------------------------------!

     pure subroutine assemble_residual_interface(this, residual)

       import :: assembler, vector

       class(assembler) , intent(inout) :: this
       class(vector)    , intent(inout) :: residual

     end subroutine assemble_residual_interface

     !----------------------------------------------------------------!
     ! Interface for jacobian assembly
     !----------------------------------------------------------------!

     pure subroutine assemble_jacobian_interface(this, jacobian)

       import :: assembler, matrix

       class(assembler) , intent(inout) :: this
       class(matrix)    , intent(inout) :: jacobian

     end subroutine assemble_jacobian_interface

  end interface

contains

end module assembler_interface
