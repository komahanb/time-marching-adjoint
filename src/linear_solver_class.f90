#include "scalar.fpp"
!=====================================================================!
! Interface for linear system solvers
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module linear_solver_class

  implicit none

  private
  public :: linear_solver

  type, abstract :: linear_solver

     ! type options
     private
     
     ! attributes
     integer :: mat_size

   contains

     ! defined procedures
     procedure :: get_mat_size
     procedure :: set_mat_size

     ! deferred procedures

  end  type linear_solver

contains

  !===================================================================!
  ! Sets the matrix size
  !===================================================================!
  
  subroutine set_mat_size(this, mat_size) 

    class(linear_solver) :: this
    type(integer) :: mat_size

    this % mat_size = mat_size

  end subroutine set_mat_size

  !===================================================================!
  ! Return the size of the linear system                                                                   !
  !===================================================================!
  
  function get_mat_size(this)

    class(linear_solver) :: this
    type(integer) :: get_mat_size

    get_mat_size = this % mat_size

  end function get_mat_size

end module linear_solver_class

!======================================================================!
! Lapack solver type
!======================================================================!

module lapack_class

  use linear_solver_class, only: linear_solver

  implicit none

  private
  public :: lapack

  type, extends(linear_solver) :: lapack
     
  end type lapack
  
  ! interface for constructing lapack
  interface lapack
     procedure constructor
  end interface lapack

  contains

    !=================================================================!
    ! Constructor for lapack solver type
    !=================================================================!

    function constructor(mat_size) result(this)

      type(lapack)  :: this
      type(integer) :: mat_size
      
      call this % set_mat_size(mat_size)

    end function constructor

  end module lapack_class
