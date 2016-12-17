#include "scalar.fpp"

!=====================================================================!
! Module implementing the Laplacian operator for heat equation
!=====================================================================!

module heat_laplacian_operator

  ! petsc modules

  use petscsys
  use petscmat

  implicit none

  private
  public  :: init_laplacian, build_laplacian, destroy_laplacian

  type(integer) :: nx   ! maximum number of x cells
  type(integer) :: ny   ! maximum number of y cells
  type(integer) :: nz   ! maximum number of z cells
  type(integer) :: ng   ! maximum number of groups
  type(integer) :: ierr ! petsc error code

  type :: laplacian_operator

     type(Mat)            :: L        ! petsc matrix for heat laplacian operator
     type(integer)              :: n        ! dimensions of matrix
     type(integer)              :: nnz      ! max number of nonzeros
     type(integer)              :: localn   ! local size on proc
     type(integer), allocatable :: d_nnz(:) ! vector of diagonal preallocation
     type(integer), allocatable :: o_nnz(:) ! vector of off-diagonal preallocation

  end type laplacian_operator

contains

  !-------------------------------------------------------------------!
  ! Perform initialization tasks for the laplacian operator           !
  !-------------------------------------------------------------------!
  
  subroutine init_laplacian(this)

    type(laplacian_operator) :: this

    ! get indices

    ! get preallocation

    ! setup lapalcian operator
    !call MatCreateMPIAIJ(PETSC_COMM_WORLD,this%localn,this%localn,PETSC_DECIDE,&
    !     & PETSC_DECIDE,PETSC_NULL_INTEGER,this%d_nnz,PETSC_NULL_INTEGER,this%o_nnz, &
    !     & this%L,ierr)
    call MatCreate(PETSC_COMM_WORLD, this % L, ierr)
    call MatSetOption(this % L, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE, ierr)
    call MatSetOption(this % L, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE, ierr)

  end subroutine init_laplacian

  !-------------------------------------------------------------------!
  ! Build the laplacian operator
  !-------------------------------------------------------------------!

  subroutine build_laplacian(this)

    type(laplacian_operator) :: this

    type(integer) :: i          ! iteration counter for x
    type(integer) :: j          ! iteration counter for y
    type(integer) :: k          ! iteration counter for z
    type(integer) :: ierr       ! Petsc error code
    type(integer) :: row_start  ! the first local row on the processor
    type(integer) :: row_finish ! the last local row on the processor
    type(integer) :: irow       ! iteration counter over row
    type(integer) :: icol       ! column index of the inserted value

    type(scalar)  :: val        ! temporary variable for nfissxs

    ! get row bounds for this processor
    call MatGetOwnershipRange(this % L, row_start, row_finish, ierr)

    ! begin iteration loops
    ROWS: do irow = row_start, row_finish-1

       ! add the diagonal entry
       icol = irow
       val  = 2.0d0
       call MatSetValue(this % L, irow, icol, val, INSERT_VALUES, ierr)

       ! add the lower diagonal entry
       icol = irow - 1
       val  = -1.0d0
       call MatSetValue(this % L, irow, icol, val, INSERT_VALUES, ierr)

       ! add the upper diagonal entry
       icol = irow - 1
       val  = -1.0d0
       call MatSetValue(this % L, irow, icol, val, INSERT_VALUES, ierr)

    end do ROWS

    ! assemble matrix 
    call MatAssemblyBegin(this % L, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(this % L, MAT_FINAL_ASSEMBLY, ierr)

    ! print out operator to file
    call print_operator(this)

  end subroutine build_laplacian

  !-------------------------------------------------------------------!
  ! Destroy laplacian operator
  !-------------------------------------------------------------------!

  subroutine destroy_laplacian(this)

    type(laplacian_operator) :: this
    
  end subroutine destroy_laplacian

  !-------------------------------------------------------------------!
  ! Determine the dimensions and number of nonzero entries in the mat
  !-------------------------------------------------------------------!

  subroutine get_indices(this)
    
    type(laplacian_operator) :: this
    
    nx = 10
    ng = 1

    ! number of nonzeros
    this % nnz = 3*nx
    
    ! dimensions
    this % n = nx
    
  end subroutine get_indices

  !-------------------------------------------------------------------!
  ! Allocate the size of the matrix
  !-------------------------------------------------------------------!

  subroutine preallocate_matrix(this)

    type(laplacian_operator) :: this

    type(integer) :: rank          ! rank of processor
    type(integer) :: sizen         ! number of procs
    type(integer) :: i             ! iteration counter for x
    type(integer) :: j             ! iteration counter for y
    type(integer) :: k             ! iteration counter for z
    type(integer) :: g             ! iteration counter for groups
    type(integer) :: h             ! energy group when doing scattering
    type(integer) :: n             ! the extent of the matrix
    type(integer) :: irow          ! row counter
    type(integer) :: row_start     ! index of local starting row
    type(integer) :: row_end       ! index of local final row
    type(integer) :: hmat_idx      ! index in matrix for energy group h

    ! get rank and max rank of procs
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, sizen, ierr)
    
    ! get local problem size
    n = this % n

    ! determine local size, divide evenly between all other procs
    this % localn = n/(sizen)

    ! add 1 more if proc id is less than mod
    if (rank < mod(n, sizen)) this % localn = this % localn + 1

    ! determine local starting row
    row_start = 0
    if (rank < mod(n,sizen)) then
       row_start = rank*(n/sizen+1)
    else
       row_start = min(mod(n,sizen)*(n/sizen+1)+(rank - mod(n,sizen))*(n/sizen),n)
    end if

    ! determine local final row
    row_end = row_start + this % localn - 1

    ! allocate counters
    if (.not. allocated(this%d_nnz)) allocate(this%d_nnz(row_start:row_end))
    if (.not. allocated(this%o_nnz)) allocate(this%o_nnz(row_start:row_end))
    this % d_nnz = 0
    this % o_nnz = 0
    
!!$    ! begin loop around local rows
!!$    ROWS: do irow = row_start,row_end
!!$
!!$       ! initialize counters 
!!$       this % d_nnz(irow) = 1 ! already add in matrix diagonal
!!$       this % o_nnz(irow) = 0
!!$
!!$       ! get location indices
!!$       call matrix_to_indices(irow, g, i, j, k)
!!$
!!$       ! begin loop over off diagonal in-scattering
!!$       NFISS: do h = 1, ng
!!$
!!$          ! cycle though if h=g
!!$          if (h == g) then
!!$             cycle
!!$          end if
!!$
!!$          ! get neighbor matrix index
!!$          call indices_to_matrix(h, i, j, k, hmat_idx)
!!$
!!$          ! record nonzero
!!$          if (((hmat_idx-1) >= row_start) .and.                        &
!!$               &   ((hmat_idx-1) <= row_end)) then
!!$             this % d_nnz(irow) = this % d_nnz(irow) + 1
!!$          else
!!$             this % o_nnz(irow) = this % o_nnz(irow) + 1
!!$          end if
!!$
!!$       end do NFISS
!!$
!!$    end do ROWS

  end subroutine preallocate_matrix

  !-------------------------------------------------------------------!
  ! indices_to_matrix takes (x,y,z,g) indices and computes location in
  ! matrix
  !-------------------------------------------------------------------!
 
  subroutine indices_to_matrix(g, i, j, k, matidx)
    
    type(integer) :: matidx         ! the index location in matrix
    type(integer) :: i               ! current x index
    type(integer) :: j               ! current y index
    type(integer) :: k               ! current z index
    type(integer) :: g               ! current group index
    
    ! compute index
    matidx = g + ng*(i - 1) + ng*nx*(j - 1) + ng*nx*ny*(k - 1)
    
  end subroutine indices_to_matrix

  !-------------------------------------------------------------------!
  ! MATRIX_TO_INDICES 
  !-------------------------------------------------------------------!
  
  subroutine matrix_to_indices(irow, g, i, j, k)

    type(integer) :: i                    ! iteration counter for x
    type(integer) :: j                    ! iteration counter for y
    type(integer) :: k                    ! iteration counter for z
    type(integer) :: g                    ! iteration counter for groups
    type(integer) :: irow                 ! iteration counter over row (0 reference)

    ! compute indices
    g = mod(irow,ng) + 1 
    i = mod(irow,ng*nx)/ng + 1
    j = mod(irow,ng*nx*ny)/(ng*nx)+ 1
    k = mod(irow,ng*nx*ny*nz)/(ng*nx*ny) + 1 
    
  end subroutine matrix_to_indices

  !-------------------------------------------------------------------!
  ! PRINT_M_OPERATOR 
  !-------------------------------------------------------------------!

  subroutine print_operator(this)

    type(laplacian_operator) :: this

    type(PetscViewer) :: viewer

    ! write out matrix in binary file (debugging)
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,'laplacemat.bin' &
         &     ,FILE_MODE_WRITE,viewer,ierr)
    call MatView(this%L,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)

  end subroutine print_operator

end module heat_laplacian_operator
