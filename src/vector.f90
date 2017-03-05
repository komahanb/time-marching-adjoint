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

  ! Define interfaces to all abstract procedures
 
  abstract interface
     
     !================================================================!
     ! Interface to adding entryies into the vector
     !================================================================!
     
     subroutine add_entry_interface(this, idx, val)

       import vector
       
       class(vector) , intent(inout) :: this
       integer       , intent(in)    :: idx
       type(scalar)  , intent(in)    :: val

     end subroutine add_entry_interface

     type(scalar) function get_entry_interface(this, idx)

       import vector
       
       class(vector) , intent(in) :: this
       type(integer) , intent(in) :: idx

     end function get_entry_interface       

  end interface

contains
  
  !===================================================================!
  ! Returns the size of the vector
  !===================================================================!
  
  pure type(integer) function get_size(this)

    class(vector), intent(in) :: this

    get_size = this % size

  end function get_size

  !===================================================================!
  ! Sets the size of the vector
  !===================================================================!

  pure subroutine set_size(this, size)

    class(vector), intent(inout) :: this  
    type(integer), intent(in)    :: size

    this % size = size

  end subroutine set_size

end module vector_interface
