#include "scalar.fpp"
!=====================================================================!
! Module that defines a matrix type that contains dense storage of
! entries
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module dense_matrix_interface

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
  
  function constructor(row_size, col_size) result(this)

    integer            :: col_size
    integer            :: row_size
    type(dense_matrix) :: this
    
    ! Set matrix dimensions
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

  subroutine add_dense_entry(this, row, col, val)

    class(dense_matrix) :: this
    type(integer)       :: col
    type(integer)       :: row
    type(scalar)        :: val

    this % vals(row, col) = val

  end subroutine add_dense_entry

  !=================================================================!
  ! Fetch the entry corresponding to the row and column
  !=================================================================!

  type(scalar) function get_dense_entry(this, row, col) result(val)

    class(dense_matrix) :: this
    type(integer) :: row
    type(integer) :: col

    val = this % vals(row, col)
    
  end function get_dense_entry

end module dense_matrix_interface
