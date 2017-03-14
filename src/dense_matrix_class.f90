#include "scalar.fpp"
!=====================================================================!
! Module that defines a matrix type that contains dense storage of
! entries
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module dense_matrix_class

  use constants, only : WP
  use matrix_interface, only: matrix

  implicit none

  private

  public :: dense_matrix
  
  ! Specialized matrix type for dense storage
  type, extends(matrix) :: dense_matrix

     type(scalar), allocatable :: vals(:,:)

   contains

     procedure :: add_entry => add_dense_entry
     procedure :: get_entry => get_dense_entry

  end type dense_matrix
  
  ! Interfaces
  interface dense_matrix
     procedure constructor
  end interface dense_matrix

contains
  
  !===================================================================!
  ! Initializes an instance of dense matrix
  !===================================================================!
    
  pure type(dense_matrix) function constructor(row_size, col_size) &
       & result(this)

    type(integer), intent(in) :: col_size
    type(integer), intent(in) :: row_size

    call this % set_row_size(row_size)
    call this % set_col_size(col_size) 

    ! Allocate space
    allocate(this % vals(this % get_row_size(), this% get_col_size()))

    ! Zero the entries
    this % vals = 0.0_WP

  end function constructor

  !=================================================================!
  ! Adding an entry to a dense matrix
  !=================================================================!

  pure subroutine add_dense_entry(this, row, col, data)
    
    class(dense_matrix), intent(inout) :: this
    type(integer)      , intent(in)    :: col
    type(integer)      , intent(in)    :: row
    type(scalar)       , intent(in)    :: data

    this % vals(row, col) = data

  end subroutine add_dense_entry

  !=================================================================!
  ! Fetch the entry corresponding to the row and column
  !=================================================================!

  pure type(scalar) function get_dense_entry(this, row, col) result(val)

    class(dense_matrix), intent(in) :: this
    type(integer)      , intent(in) :: row
    type(integer)      , intent(in) :: col

    val = this % vals(row, col)
    
  end function get_dense_entry

end module dense_matrix_class
