program main

  ! import petsc modules

  use petscsys
  use petscvec
  use petscmat
  use petscpc
  use petscksp

  implicit none

  type(vec)     :: x, b, u
  type(mat)     :: A
  type(ksp)     :: cksp

  type(real)    :: norm


  type(integer) :: dim, i, j, ii, jj, istart, iend, its
  type(logical) :: flg
  type(integer) :: n = 100
  type(integer) :: ierr, rank, size
  type(integer) :: irow, icol

  type(PetscViewer)  :: viewer
!  type(PetscRandom)  :: rctx

  print *, "Initializing petsc"
  
  call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

!  call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
!  call MPI_Comm_size(PETSC_COMM_WORLD, size, ierr)

  ! Create matrix
  call MatCreate(PETSC_COMM_WORLD, A, ierr)
  call MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n, ierr)
  call MatSetFromOptions(A, ierr)
  call MatSetUp(A, ierr)

  call MatGetOwnershipRange(A, istart, iend, ierr)
 
  print *, istart, iend
       
  ! begin iteration loops
  ROWS: do irow = istart, iend - 1

     icol = 1

     call MatSetValue(A,irow,icol,1.0,INSERT_VALUES,ierr)

  end do ROWS

  ! assemble matrix 
  call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, ierr)
  call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, ierr)
  
  ! write out matrix in binary file (debugging)
  call PetscViewerBinaryOpen(PETSC_COMM_WORLD,'mat.bin', FILE_MODE_WRITE, viewer, ierr)
  call MatView(A,viewer,ierr)
  call PetscViewerDestroy(viewer,ierr)

!!$  ! create vectors
!!$  call VecCreate(PETSC_COMM_WORLD, x, ierr)
!!$  call VecSetSizes(x, PETSC_DECIDE, n, ierr)
!!$  !call VecSetFromOptions(x, ierr)
!!$
!!$  call VecCreate(PETSC_COMM_WORLD, u, ierr)
!!$  call VecSetSizes(u, PETSC_DECIDE, n, ierr)
!!$  !call VecSetFromOptions(u, ierr)
!!$
!!$  call VecCreate(PETSC_COMM_WORLD, b, ierr)
!!$  call VecSetSizes(b, PETSC_DECIDE, n, ierr)
!!$  !call VecSetFromOptions(b, ierr)

  !-------------------------------------------------------------------!
  ! Create the linear solver
  !-------------------------------------------------------------------!

  call KSPCreate(PETSC_COMM_WORLD, cksp, ierr)
  call KSPSetOperators(cksp,A,A,DIFFERENT_NONZERO_PATTERN,ierr)

  ! The matrix and preconditioner are the same

!  call KSPGetPC(cksp, cpc, ierr)
!  ptype = PCJACOBI
!  call PCSetType(pc,ptype,ierr)
!      tol = 1.e-7
!      call KSPSetTolerances(ksp,tol,PETSC_DEFAULT_DOUBLE_PRECISION,
!     &     PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_INTEGER,ierr)
! call KSPSetFromOptions(cksp,ierr)
 
 !-------------------------------------------------------------------!
 ! Solve the linear system
 !-------------------------------------------------------------------!
 
 call KSPSolve(cksp, b, x, ierr)

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

! call KSPDestroy(cksp,ierr)
! call VecDestroy(u,ierr)
! call VecDestroy(x,ierr)
! call VecDestroy(b,ierr)
! call MatDestroy(A,ierr)

 print *, "Finalizing petsc done"
 call PetscFinalize(ierr)

end program main
