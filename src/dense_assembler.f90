#include "scalar.fpp"
!=====================================================================!
! Module that implements a derived type called 'dense_assembler' which
! produces a dense matrix and vector that goes into linear and
! nonlinear solution routines.
! 
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module dense_assembler_class

  use assembler_interface    , only : assembler
  use matrix_interface  , only : matrix
  use vector_interface  , only : vector
  use dense_matrix_interface , only : dense_matrix
  use dense_vector_interface , only : dense_vector
  
  implicit none
  
  type, extends(assembler) :: dense_assembler

   contains

     ! Implement deferred procedures from abstract type
     procedure :: assemble_residual
     procedure :: assemble_jacobian

  end type dense_assembler

contains

  !===================================================================!
  ! Implementation for residual assembly
  !===================================================================!
  
  pure subroutine assemble_residual(this, residual)

    class(dense_assembler) , intent(inout) :: this
    class(vector)          , intent(inout) :: residual

  end subroutine assemble_residual

  !===================================================================!
  ! Implementation for jacobian assembly
  !===================================================================!

  pure subroutine assemble_jacobian(this, jacobian)

    class(dense_assembler) , intent(inout) :: this
    class(matrix)          , intent(inout) :: jacobian

  end subroutine assemble_jacobian
  
end module dense_assembler_class
