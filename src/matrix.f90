#include "scalar.fpp"
!=====================================================================!
! Module that contains matrix type definition and procedures to add
! and get entries from the matrix. This class also defines the
! interface to add entries into the matrix and get entries from the
! matrix.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module matrix_interface

  implicit none

  private
  public :: matrix

  type, abstract :: matrix
     
     private

     ! the number of rows in a matrix
     type(integer) :: row_size

     ! the number of columns in a matrix
     type(integer) :: col_size 

   contains

     ! defined procedures
     procedure :: get_col_size
     procedure :: get_row_size
     procedure :: set_col_size
     procedure :: set_row_size
     
     ! deferred procedures
     procedure(add_entry_interface), deferred :: add_entry
     procedure(get_entry_interface), deferred :: get_entry
     
  end type matrix

  ! Define interfaces to all abstract procedures
  abstract interface

     !----------------------------------------------------------------!     
     ! Adding entry into a matrix  
     !----------------------------------------------------------------!

     pure subroutine add_entry_interface(this, row, col, data)

       import :: matrix

       class(matrix), intent(inout) :: this
       type(integer), intent(in)    :: row
       type(integer), intent(in)    :: col
       type(scalar) , intent(in)    :: data

     end subroutine add_entry_interface

     !----------------------------------------------------------------!
     ! Getting a scalar entry from a matrix
     !----------------------------------------------------------------!
     
     pure type(scalar) function get_entry_interface(this, row, col)

       import :: matrix

       class(matrix), intent(in) :: this
       type(integer), intent(in) :: row
       type(integer), intent(in) :: col

     end function get_entry_interface

  end interface

contains
  
  !===================================================================!
  ! Returns the column size
  !===================================================================!

  pure type(integer) function get_col_size(this)

    class(matrix), intent(in) :: this
    
    get_col_size = this % col_size

  end function get_col_size

  !===================================================================!
  ! Returns the row size
  !===================================================================!

  pure type(integer) function get_row_size(this)

    class(matrix), intent(in) :: this

    get_row_size = this % row_size

  end function get_row_size

  !===================================================================!
  ! Sets the row size
  !===================================================================!

  pure subroutine set_row_size(this, row)

    class(matrix), intent(inout) :: this
    type(integer), intent(in)    :: row

    this % row_size = row
    
  end subroutine set_row_size
  
  !===================================================================!
  ! sets the column size
  !===================================================================!
  
  pure subroutine set_col_size(this, col)
    
    class(matrix), intent(inout) :: this
    type(integer), intent(in)    :: col

    this % col_size = col

  end subroutine set_col_size
  
end module matrix_interface
