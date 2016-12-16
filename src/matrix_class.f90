#include "scalar.fpp"

!=====================================================================!
! Module that contains matrix type definition and procedures to add
! and get entries from the matrix. This class also defines the
! interface to add elements into the matrix.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module matrix_class

  implicit none

  private
  public :: matrix

  type, abstract :: matrix
     
     ! type options
     private

     ! attributes
     integer :: row_size ! the number of rows in a matrix
     integer :: col_size ! the number of columns in a matrix

   contains

     ! defined procedures
     procedure :: get_col_size
     procedure :: get_row_size
     procedure :: set_col_size
     procedure :: set_row_size
     
     ! deferred procedures
     procedure(add_element_interface), deferred :: add_element
     
  end type matrix

  ! define interfaces to all abstract procedures
 
  abstract interface
     
     !================================================================!
     ! interface to adding elements into the matrix
     !================================================================!
     
     subroutine add_element_interface(this, row, col, val)

       import matrix

       class(matrix) :: this
       type(integer) :: row
       type(integer) :: col
       type(scalar)  :: val

     end subroutine add_element_interface

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
  
end module matrix_class
