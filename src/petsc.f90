#include "scalar.fpp"

program main

  ! import petsc modules

  use petscsys
  use petscvec
  use petscmat
  use petscpc
  use petscksp

  implicit none

  type(integer)     :: dim
  type(integer)     :: n = 100
  type(integer)     :: i, j, ii, jj 
  type(integer)     :: istart, iend
  type(integer)     :: iters
  type(integer)     :: ierr, rank, size
  type(integer)     :: irow, icol

  type(logical)     :: flag
  type(real)        :: norm

  type(vec)         :: x, b, ax
  type(mat)         :: A
  type(ksp)         :: cksp
  type(pc)          :: cpc

  type(PetscViewer) :: viewer
  type(scalar)      :: tmp

  print *, "Initializing petsc"
  
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

  ! Create matrix
  call MatCreate(PETSC_COMM_WORLD, A, ierr)
  call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n, ierr)
  !call MatSetFromOptions(A, ierr)
  call MatSetUp(A, ierr)

  ! Get the ownership range of the matrix on each processor
  call MatGetOwnershipRange(A, istart, iend, ierr)

  ! print *, "ownership of rank", rank, " between ", istart, iend
  
  ! set entries into the matrix

  ! Create a central difference operator matrix
  rows: do irow = istart, iend - 1
     cols: do icol = istart, iend - 1
        if (irow .eq. icol) then
           tmp = 2.0d0  ! d
        else if (irow .eq. icol - 1) then
           tmp = -1.0d0 ! -r
        else if (icol .eq. icol + 1) then
           tmp = -1.0d0 ! -r
        end if
        call MatSetValue(A, irow, icol, tmp, INSERT_VALUES, ierr)
     end do cols
  end do rows

  ! assemble matrix ( once the matrix is assembled, the matrix can
  ! not be changed )
  
  call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)

  ! Optionally write and read matrix in binary file (debugging)
  ! call PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'mat.bin', FILE_MODE_WRITE, viewer, ierr)
  ! call PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'mat.bin', FILE_MODE_READ, viewer, ierr)
  ! call MatView(A, viewer, ierr)
  ! call PetscViewerDestroy(viewer,ierr)  

  ! create vectors Create parallel vectors.  - Here, the parallel
!  partitioning of the vector is determined by PETSc at runtime.  We
!  could also specify the local dimensions if desired -- or use the
!  more general routine VecCreate().  - When solving a linear system,
!  the vectors and matrices MUST be partitioned accordingly.  PETSc
!  automatically generates appropriately partitioned matrices and
!  vectors when MatCreate() and VecCreate() are used with the same
!  communicator.  - Note: We form 1 vector from scratch and then
!  duplicate as needed.

  call VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, n, ax, ierr)
  call VecSetFromOptions(ax, ierr)
  call VecDuplicate(ax, b, ierr)
  call VecDuplicate(b, x, ierr)

!!$
!!$  ! rhs vector
!!$  call VecCreate(PETSC_COMM_WORLD, b, ierr)
!!$  call VecSetSizes(b, PETSC_DECIDE, n, ierr)
!!$  call VecSetFromOptions(b, ierr)
!!$  
!!$  ! solution vector
!!$  call VecCreate(PETSC_COMM_WORLD, x, ierr)
!!$  call VecSetSizes(x, PETSC_DECIDE, n, ierr)
!!$  call VecSetFromOptions(x, ierr)
!!$
!!$  ! temporary vector to store Ax
!!$  call VecCreate(PETSC_COMM_WORLD, ax, ierr)
!!$  call VecSetSizes(ax, PETSC_DECIDE, n, ierr)
!!$  call VecSetFromOptions(ax, ierr)

!  call VecSet(b, 1.0d0, ierr)
!  call MatMult(A, x, b)

  !-------------------------------------------------------------------!
  ! Create the linear solver
  !-------------------------------------------------------------------!
 !  Set exact solution; then compute right-hand-side vector.
!  By default we use an exact solution of a vector with all
!  elements of 1.0;  Alternatively, using the runtime option
!  -random_sol forms a solution vector with random components.

!      call PetscOptionsHasName(PETSC_NULL_CHARACTER,                    &
!     &             "-random_exact_sol",flg,ierr)
!      if (flg) then
 !        call PetscRandomCreate(PETSC_COMM_WORLD,rctx,ierr)
 !        call PetscRandomSetFromOptions(rctx,ierr)
 !        call VecSetRandom(u,rctx,ierr)
 !        call PetscRandomDestroy(rctx,ierr)
 !     else
  call VecSet(ax, 1.0d0, ierr)
  !     endif
  call MatMult(A, ax, b, ierr)

!  View the exact solution vector if desired

!      call PetscOptionsHasName(PETSC_NULL_CHARACTER,                    &
!     &             "-view_exact_sol",flg,ierr)
 !     if (flg) then
!  call VecView(ax, PETSC_VIEWER_STDOUT_WORLD, ierr)

  !    endif

  call KSPCreate(PETSC_COMM_WORLD, cksp, ierr)

  ! currently the matrix itself serves as the preconditioner
  !  call KSPSetOperators(cksp, A, A, ierr) 
  call KSPSetOperators(cksp, A, A, DIFFERENT_NONZERO_PATTERN, ierr)

  ! The matrix and preconditioner are the same
!  call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-my_ksp_monitor',  &
!       &                    flg,ierr)
!  if (flg) then
!  call KSPMonitorSet(sksp, MyKSPMonitor, PETSC_NULL_OBJECT,          &
!       &                     PETSC_NULL_FUNCTION, ierr)
  !  endif

 !-------------------------------------------------------------------!
 ! Solve the linear system
 !-------------------------------------------------------------------!
 
  call KSPSolve(cksp, b, x, ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                     Check solution and clean up
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  !  Check the error
  call VecAXPY(x, -1.0d0, b, ierr)
  call VecNorm(x, NORM_2, norm, ierr)
  call KSPGetIterationNumber(cksp, iters, ierr)

  if (norm .gt. 1.e-12) then
     print *, norm, iters
  else
     print *, norm, iters
  endif

  call VecView(ax, PETSC_VIEWER_STDOUT_WORLD, ierr)
  call VecView(b, PETSC_VIEWER_STDOUT_WORLD, ierr)

!!$ ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$ !                     Check solution and clean up
!!$ ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$
!!$ !  Check the error
!!$ call VecAXPY(x,-1.0d0,u,ierr)
!!$ call VecNorm(x,NORM_2,norm,ierr)
!!$ call KSPGetIterationNumber(cksp,its,ierr)
!!$ if (rank .eq. 0) then
!!$    if (norm .gt. 1.e-12) then
!!$       write(6,100) norm,its
!!$    else
!!$       write(6,110) its
!!$    endif
!!$ endif
!!$100 format('Norm of error ',e11.4,' iterations ',i5)
!!$110 format('Norm of error < 1.e-12,iterations ',i5)

 !  Free work space.  All PETSc objects should be destroyed when they
 !  are no longer needed.

  print *, "Destroy KSP"
! call KSPDestroy(cksp, ierr)

 print *, "Destroy vectors"
! call VecDestroy(ax, ierr)
 print *, "Destroy vectors"
! call VecDestroy(x, ierr)
 print *, "Destroy vectors"
! call VecDestroy(b, ierr)

! call MatDestroy(A, ierr)

 print *, "Finalizing petsc done"
 call PetscFinalize(ierr)
 print *, "Execution complete"
end program main
