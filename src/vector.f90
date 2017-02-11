#include "scalar.fpp"

!=====================================================================!
! Module that contains vector type definition and procedures to add
! and get entries from the vector. This class also defines the
! interface to add elements into the vector.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module vector_interface

  implicit none

  private
  public :: vector

  type, abstract :: vector
     
     ! type options
     private

     ! attributes
     integer :: size ! the number of columns in a vector

   contains

     ! defined procedures
     procedure :: get_size
     procedure :: set_size
     
     ! deferred procedures
     procedure(add_entry_interface), deferred :: add_entry
     procedure(get_entry_interface), deferred :: get_entry

  end type vector

  ! define interfaces to all abstract procedures
 
  abstract interface
     
     !================================================================!
     ! interface to adding entryies into the vector
     !================================================================!
     
     subroutine add_entry_interface(this, idx, val)

       import vector
       
       class(vector) :: this
       type(integer) :: idx
       type(scalar)  :: val

     end subroutine add_entry_interface

     type(scalar) function get_entry_interface(this, idx)

       import vector
       
       class(vector) :: this
       type(integer) :: idx

     end function get_entry_interface
       

  end interface

contains
  
  !===================================================================!
  ! returns the size
  !===================================================================!
  
  function get_size(this)

    class(vector) :: this
    type(integer) :: get_size

    get_size = this % size

  end function get_size

  !===================================================================!
  ! sets the size
  !===================================================================!

  subroutine set_size(this, size)

    class(vector) :: this
    type(integer) :: size

    this % size = size

  end subroutine set_size

end module vector_interface
