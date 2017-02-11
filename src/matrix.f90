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

     subroutine add_entry_interface(this, row, col, val)

       import :: matrix

       class(matrix) :: this
       type(integer) :: row
       type(integer) :: col
       type(scalar)  :: val

     end subroutine add_entry_interface

     !----------------------------------------------------------------!
     ! Getting a scalar entry from a matrix
     !----------------------------------------------------------------!
     
     type(scalar) function get_entry_interface(this, row, col)

       import :: matrix

       class(matrix) :: this
       type(integer) :: row
       type(integer) :: col

     end function get_entry_interface

  end interface

contains
  
  !===================================================================!
  ! returns the column size
  !===================================================================!

  function get_col_size(this)

    class(matrix) :: this
    type(integer) :: get_col_size
    
    get_col_size = this % col_size

  end function get_col_size

  !===================================================================!
  ! returns the row size
  !===================================================================!

  function get_row_size(this)

    class(matrix) :: this
    type(integer) :: get_row_size

    get_row_size = this % row_size

  end function get_row_size

  !===================================================================!
  ! sets the row size
  !===================================================================!

  subroutine set_row_size(this, row)

    class(matrix) :: this
    type(integer) :: row

    this % row_size = row
    
  end subroutine set_row_size
  
  !===================================================================!
  ! sets the column size
  !===================================================================!
  
  subroutine set_col_size(this, col)
    
    class(matrix) :: this
    type(integer) :: col

    this % col_size = col

  end subroutine set_col_size
  
end module matrix_interface
