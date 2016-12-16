module sparse_matrix_class

  ! module references
  use matrix_class,  only: matrix

  ! module options
  implicit none
  private
  public :: sparse_matrix

  ! module definitions
  type, extends(matrix) :: sparse_matrix

    ! type options
    private

    ! type attributes
    integer, allocatable :: cols(:)
    integer, allocatable :: rows(:)
    real(8), allocatable :: vals(:)

    ! type procedure definitions
    contains

      procedure :: add_element => add_sparse_element

  end type sparse_matrix

  ! interfaces
  interface sparse_matrix
    procedure constructor ! add constructor to petsc_matrix generic interface
  end interface sparse_matrix

  ! procedures
  contains

!===============================================================================
! initializes an instance of petsc_matrix
!===============================================================================

    function constructor(row_size, col_size) result(this)

      ! arguments
      integer :: col_size
      integer :: row_size
      type(sparse_matrix) :: this

      ! begin execution
      call this % set_row_size(row_size)
      call this % set_col_size(col_size) 

      print *, 'I AM CONSTRUCTING'

    end function constructor

!===============================================================================
! adding an element to a petsc matrix
!===============================================================================

    subroutine add_sparse_element(this,row,col,val)

      ! arguments
      class(sparse_matrix) :: this
      integer :: col
      integer :: row
      real(8) :: val

      ! print that we did this
      print *, 'Added values to sparse matrix!'

    end subroutine add_sparse_element

end module sparse_matrix_class
