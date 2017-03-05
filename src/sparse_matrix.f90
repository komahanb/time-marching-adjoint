#include "scalar.fpp"
!=====================================================================!
! A sparse matrix type for sparse linear algebra
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module sparse_matrix_interface

  use matrix_interface, only: matrix

  implicit none

  private
 
  public :: sparse_matrix
  
  ! Specialized matrix type for sparse storage
  type, extends(matrix) :: sparse_matrix

     type(integer), allocatable :: indices(:)
     type(integer), allocatable :: indptr(:)
     type(scalar) , allocatable :: data(:)

   contains

      procedure :: add_entry => add_sparse_entry
      procedure :: get_entry => get_sparse_entry

  end type sparse_matrix

  ! Interfaces
  interface sparse_matrix
     procedure constructor
  end interface sparse_matrix

  contains

    !=================================================================!
    ! Initializes an instance of sparse matrix
    !=================================================================!
    
    pure type(sparse_matrix) function constructor(row_size, col_size) &
         & result(this)

      type(integer), intent(in) :: col_size
      type(integer), intent(in) :: row_size

      call this % set_row_size(row_size)
      call this % set_col_size(col_size) 
      
      ! stop "SPARSE_MATRIX: Unimplemented"

    end function constructor

    !=================================================================!
    ! Adding an entry to a sparse matrix
    !=================================================================!

    pure subroutine add_sparse_entry(this, row, col, data)

      class(sparse_matrix), intent(inout) :: this
      type(integer)       , intent(in)    :: col
      type(integer)       , intent(in)    :: row
      type(scalar)        , intent(in)    :: data

      !print *, 'Added values to sparse matrix!'

    end subroutine add_sparse_entry

    !=================================================================!
    ! Fetch the entry corresponding to the row and column
    !=================================================================!
    
    pure type(scalar) function get_sparse_entry(this, row, col)

      class(sparse_matrix) , intent(in) :: this
      type(integer)        , intent(in) :: row
      type(integer)        , intent(in) :: col

    end function get_sparse_entry
    
end module sparse_matrix_interface
