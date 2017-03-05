#include "scalar.fpp"

!=====================================================================!
! Dense vector implementation of vector interface. 
! 
! TODO: Unlimited polymorphic entries to support any datatype and
! block vectors. Currently the vectors are assumed to be of
! type(scalar)
! 
! Author: Komahan Boopathy (komahan@gatech.edu)
!
!=====================================================================!

module dense_vector_interface

  use constants        , only : WP
  use vector_interface , only : vector

  implicit none
  
  private
  
  public :: dense_vector

  ! Specialized vector type
  type, extends(vector) :: dense_vector
     
     type(scalar), allocatable :: vals(:)

     contains
       
       procedure :: add_entry => add_dense_entry
       procedure :: get_entry => get_dense_entry

  end type dense_vector

  ! Interface
  interface dense_vector
     procedure constructor
  end interface dense_vector
  
contains
  
  !===================================================================!
  ! Constructor for dense matrix
  !===================================================================!

  pure type(dense_vector) function constructor (size) result (this)

    integer, intent(in) :: size

    ! Set vector dimensions
    call this % set_size(size)

    ! Allocate space
    allocate(this % vals(this % get_size()))

    ! Zero the entries
    this % vals = 0.0_WP

  end function constructor

  !===================================================================!
  ! Add an entry into the vector's specified index
  !===================================================================!

  pure subroutine add_dense_entry(this, idx, val)

    class(dense_vector), intent(inout) :: this
    integer            , intent(in)    :: idx
    type(scalar)       , intent(in)    :: val

    this % vals(idx) = val

  end subroutine add_dense_entry

  !===================================================================!
  ! Get entry of the vector at the supplied index
  !===================================================================!
  
  pure type(scalar) function get_dense_entry(this, idx) result(val)

    class(dense_vector) :: this
    type(integer)       :: idx

    val = this % vals(idx)

  end function get_dense_entry

end module dense_vector_interface
