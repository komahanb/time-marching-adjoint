!=====================================================================!
! Module that contains handy linear algebra operations using LAPACK.
!=====================================================================!

module linear_algebra
  
  ! import dependencies
  use iso_fortran_env, only : dp => REAL64
  
  use lapack, only: dsyevd, dsygvd, ilaenv, zgetri, zgetrf, zheevd, &
                    dgeev , zgeev , zhegvd, dgesv , zgesv , dgetrf, &
                    dgetri, dgelsy, zgelsy, dgesvd, zgesvd
  
  ! disable implicit datatypes
  implicit none

  ! restrict access to all functions
  private

  ! expose only a required functions
  public :: eig, eigvals, eigh              ! eigenvals and vec
  public :: inv, solve, lstsq, svdvals, svd ! linear system
  public :: det, eye, diag, trace           ! matrix utils
  public :: assert, stop_error              ! misc
  
  ! Module parameters  
  complex(dp), parameter :: i_ = (0, 1)
  
  !===================================================================!
  ! Eigen interfaces
  !===================================================================!

  ! eigenvalue/-vector problem for general matrices:
  interface eig
     module procedure deig
     module procedure zeig
  end interface eig

  ! eigenvalue/-vector problem for real symmetric/complex hermitian
  ! matrices:
  interface eigh
     module procedure deigh_generalized
     module procedure deigh_simple
     module procedure zeigh_generalized
     module procedure zeigh_simple
  end interface eigh

  ! eigenvalues for general matrices:
  interface eigvals
     module procedure deigvals
     module procedure zeigvals
  end interface eigvals

  !===================================================================!
  ! Linear solve interfaces
  !===================================================================!

  ! matrix inversion for real/complex matrices:
  interface inv
     module procedure dinv
     module procedure zinv
  end interface inv

  ! solution to linear systems of equation with real/complex
  ! coefficients:
  interface solve
     module procedure dsolve
     module procedure zsolve
  end interface solve

  ! least square solutions the real/complex systems of equations of
  ! possibly non-square shape:
  interface lstsq
     module procedure dlstsq
     module procedure zlstsq
  end interface lstsq

  ! singular values of real/complex matrices:
  interface svdvals
     module procedure dsvdvals
     module procedure zsvdvals
  end interface svdvals

  ! singular value decomposition of real/complex matrices:
  interface svd
     module procedure dsvd
     module procedure zsvd
  end interface svd

  !===================================================================!
  ! Utility interfaces
  !===================================================================!
  
  ! determinants of real/complex square matrices:
  interface det
     module procedure ddet
     module procedure zdet
  end interface det

  ! construction of square matrices from the diagonal elements:
  interface diag
     module procedure ddiag
     module procedure zdiag
  end interface diag

  ! trace of real/complex matrices:
  interface trace
     module procedure dtrace
     module procedure ztrace
  end interface trace

  ! assert shape of matrices:
  interface assert_shape
     module procedure dassert_shape
     module procedure zassert_shape
  end interface assert_shape

contains

  !-------------------------------------------------------------------!
  ! Get the eigenvalues and eigenvectors of a real matrix.
  !
  ! TODO: add optional switch for left or right eigenvectors in deig()
  ! and zeig()?
  ! -------------------------------------------------------------------!
  
  subroutine deig(A, lam, c)

    real(dp), intent(in)     :: A(:,:)  ! matrix for eigenvalue compuation
    complex(dp), intent(out) :: lam(:)  ! eigenvalues: A c = lam c
    complex(dp), intent(out) :: c(:,:)  ! eigenvectors: A c = lam c; c(i,j) = ith component of jth vec.

    ! LAPACK variables for DGEEV:
    real(dp), allocatable :: At(:,:), vl(:,: ), vr(:,:)
    real(dp), allocatable :: wi(:), work(:), wr(:)
    integer               :: info, lda, ldvl, ldvr, lwork, n, i
    
    preprocess: block

      ! Determine the size of the linear system
      lda = size(A(:,1))
      n = size(A(1,:))
      call assert_shape(A, [n, n], "solve", "A")
      call assert_shape(c, [n, n], "solve", "c")
      ldvl = n
      ldvr = n
      lwork = 8*n  ! TODO: can this size be optimized? query first?
      allocate(At(lda,n), wr(n), wi(n), vl(ldvl,n), vr(ldvr,n), work(lwork))
      At = A

    end block preprocess

    solve: block

      ! Call LAPACK to find the eigen values and vectors
      call dgeev('N', 'V', n, At, lda, wr, wi, vl, ldvl, vr, ldvr, &
           work, lwork, info)

    end block solve

    postprocess: block

      if (info .ne. 0) then
         print *, "dgeev returned info = ", info
         if (info .lt. 0) then
            print *, "the", -info, "-th argument had an illegal value"
         else
            print *, "the QR algorithm failed to compute all the"
            print *, "eigenvalues, and no eigenvectors have been computed;"
            print *, "elements ", info+1, ":", n, "of WR and WI contain eigenvalues which"
            print *, "have converged."
         end if
         call stop_error('eig: dgeev error')
      end if

      lam = wr + i_*wi
      ! As DGEEV has a rather complicated way of returning the
      ! eigenvectors, it is necessary to build the complex array of
      ! eigenvectors from two real arrays:
      do concurrent (i = 1: n)
         if (wi(i) > 0.0) then  
            ! first of two conjugate eigenvalues
            c(:, i) = vr(:, i) + i_*vr(:, i+1)
         elseif (wi(i) < 0.0_dp) then
            ! second of two conjugate eigenvalues
            c(:, i) = vr(:, i-1) - i_*vr(:, i)
         else
            c(:, i) = vr(:, i)
         end if
      end do

    end block postprocess

  end subroutine deig

  !-------------------------------------------------------------------!
  ! Get the eigenvalues and eigenvectors of a complex matrix.
  !
  ! TODO: add optional switch for left or right eigenvectors in deig()
  ! and zeig()?
  !-------------------------------------------------------------------!
  
  subroutine zeig(A, lam, c)

    complex(dp), intent(in)  :: A(:, :) ! matrix to solve eigenproblem for
    complex(dp), intent(out) :: lam(:)  ! eigenvalues: A c = lam c
    complex(dp), intent(out) :: c(:,:)  ! eigenvectors: A c = lam c;
                                        ! c(i,j) = ith component of
                                        ! jth vec.

    ! LAPACK variables: 
    integer                  :: info, lda, ldvl, ldvr, lwork, n, lrwork 
    real(dp), allocatable    :: rwork(:) 
    complex(dp), allocatable :: vl(:,:), vr(:,:), work(:)
    
    preprocess: block

      ! Determine the size of matrix
      lda = size(A(:,1))
      n = size(A(1,:))
      call assert_shape(A, [n, n], "solve", "A")
      call assert_shape(c, [n, n], "solve", "c")
      ldvl = n
      ldvr = n
      lwork = 8*n  ! TODO: can this size be optimized? query first?
      lrwork = 2*n
      allocate(vl(ldvl,n), vr(ldvr,n), work(lwork), rwork(lrwork))
      c = A

    end block preprocess

    solve: block

      ! Find the eigen values
      call zgeev('N', 'V', n, c, lda, lam, vl, ldvl, vr, ldvr, work, &
           lwork, rwork, info)

    end block solve

    postprocess: block

      if (info .ne. 0) then
         print *, "zgeev returned info = ", info
         if (info .lt. 0) then
            print *, "the ",-info, "-th argument had an illegal value."
         else
            print *, "the QR algorithm failed to compute all the"
            print *, "eigenvalues, and no eigenvectors have been computed;"
            print *, "elements and ", info+1, ":", n, " of W contain eigenvalues which have"
            print *, "converged."
         end if
         call stop_error('eig: zgeev error')
      end if

      c = vr

    end block postprocess

  end subroutine zeig

  !-------------------------------------------------------------------!
  ! Returns the eigen values of the matrix
  !-------------------------------------------------------------------!
  
  function deigvals(A) result(lam)

    real(dp), intent(in)     :: A(:, :) ! matrix for eigenvalue compuation
    complex(dp), allocatable :: lam(:)  ! eigenvalues: A c = lam c
    
    ! LAPACK variables for DGEEV:
    real(dp), allocatable :: At(:,:), vl(:,: ), vr(:,:)
    real(dp), allocatable :: wi(:), work(:), wr(:)
    integer               :: info, lda, ldvl, ldvr, lwork, n


    preprocess: block

      ! Determine the size of the matrix
      lda = size(A(:,1))
      n = size(A(1,:))
      call assert_shape(A, [n, n], "solve", "A")
      ldvl = n
      ldvr = n
      lwork = 8*n  ! TODO: can this size be optimized? query first?
      allocate(At(lda,n), wr(n), wi(n), vl(ldvl,n), vr(ldvr,n), work(lwork), lam(n))
      At = A

    end block preprocess

    solve: block

      ! Determine the eigenvalues
      call dgeev('N', 'N', n, At, lda, wr, wi, vl, ldvl, vr, ldvr, &
           work, lwork, info)

    end block solve

    postprocess: block

      ! Post process
      if (info .ne. 0) then
         print *, "dgeev returned info = ", info
         if (info .lt. 0) then
            print *, "the", -info, "-th argument had an illegal value"
         else
            print *, "the QR algorithm failed to compute all the"
            print *, "eigenvalues, and no eigenvectors have been computed;"
            print *, "elements ", info+1, ":", n, "of WR and WI contain eigenvalues which"
            print *, "have converged."
         end if
         call stop_error('eigvals: dgeev error')
      end if

      lam = wr + i_*wi

    end block postprocess

  end function deigvals

  !===================================================================!
  ! Eigenvalues of a complex matrix
  !===================================================================!

  function zeigvals(A) result(lam)

    complex(dp), intent(in)  :: A(:, :)  ! matrix to solve eigenproblem for
    complex(dp), allocatable :: lam(:)  ! eigenvalues: A c = lam c

    ! LAPACK variables:
    integer                  :: info, lda, ldvl, ldvr, lwork, n, lrwork
    real(dp), allocatable    :: rwork(:)
    complex(dp), allocatable :: At(:,:), vl(:,:), vr(:,:), work(:)

    preprocess: block

      ! Determine the size of the matrix
      lda = size(A(:,1))
      n = size(A(1,:))
      call assert_shape(A, [n, n], "solve", "A")
      ldvl = n
      ldvr = n
      lwork = 8*n  ! TODO: can this size be optimized? query first?
      lrwork = 2*n
      allocate(At(lda,n), &
           & vl(ldvl,n), vr(ldvr,n), &
           & work(lwork), rwork(lrwork), lam(n))
      At = A

    end block preprocess

    solve: block

      ! Determine the eigen values
      call zgeev('N', 'N', n, At, lda, lam, vl, ldvl, vr, ldvr, work, &
           lwork, rwork, info)

    end block solve

    postprocess: block

      ! Post process
      if (info .ne. 0) then
         print *, "zgeev returned info = ", info
         if (info .lt. 0) then
            print *, "the ",-info, "-th argument had an illegal value."
         else
            print *, "the QR algorithm failed to compute all the"
            print *, "eigenvalues, and no eigenvectors have been computed;"
            print *, "elements and ", info+1, ":", n, " of W contain eigenvalues which have"
            print *, "converged."
         end if
         call stop_error('eig: zgeev error')
      end if

    end block postprocess

  end function zeigvals

  !-------------------------------------------------------------------!
  ! Solves generalized eigen value problem for all eigenvalues and
  ! eigenvectors. Am must be symmetric and Bm symmetric positive
  ! definite. Only the lower triangular part of Am and Bm is used in
  ! computations.
  !-------------------------------------------------------------------!
  
  subroutine deigh_generalized(Am, Bm, lam, c)

    real(dp), intent(in)  :: Am(:,:)  ! LHS matrix: Am c = lam Bm c
    real(dp), intent(in)  :: Bm(:,:)  ! RHS matrix: Am c = lam Bm c
    real(dp), intent(out) :: lam(:)   ! eigenvalues: Am c = lam Bm c
    real(dp), intent(out) :: c(:,:)   ! eigenvectors: Am c = lam Bm c;
                                      ! c(i,j) = ith component of jth
                                      ! vec.
    integer :: n
    
    ! LLAPACK variables
    integer               :: lwork, liwork, info
    integer, allocatable  :: iwork(:)
    real(dp), allocatable :: Bmt(:,:), work(:)

    preprocess: block
      
      ! Matrix size
      n = size(Am,1)
      call assert_shape(Am, [n, n], "eigh", "Am")
      call assert_shape(Bm, [n, n], "eigh", "B")
      call assert_shape(c, [n, n], "eigh", "c")

      ! Determine the work size arrays
      lwork = 1 + 6*n + 2*n**2
      liwork = 3 + 5*n
      allocate(Bmt(n,n), work(lwork), iwork(liwork))
      c = Am; Bmt = Bm  ! Bmt temporaries overwritten by dsygvd

    end block preprocess

    solve: block
      
      call dsygvd(1,'V','L',n,c,n,Bmt,n,lam,work,lwork,iwork,liwork,info)

    end block solve

    postprocess: block
      
      ! Post process
      if (info .ne. 0) then
         print *, "dsygvd returned info =", info
         if (info .lt. 0) then
            print *, "the", -info, "-th argument had an illegal value"
         else if (info .le. n) then
            print *, "the algorithm failed to compute an eigenvalue while working"
            print *, "on the submatrix lying in rows and columns", 1.0_dp*info/(n+1)
            print *, "through", mod(info, n+1)
         else
            print *, "The leading minor of order ", info-n, &
                 "of B is not positive definite. The factorization of B could ", &
                 "not be completed and no eigenvalues or eigenvectors were computed."
         end if
         call stop_error('eigh: dsygvd error')
      end if

    end block postprocess

  end subroutine deigh_generalized

  !-------------------------------------------------------------------!
  ! Solves eigenvalue problem for all eigenvalues and eigenvectors. Am
  ! must by symmetric. Only the lower triangular part of Am is used.
  !-------------------------------------------------------------------!
  
  subroutine deigh_simple(Am, lam, c)

    real(dp), intent(in)  :: Am(:,:)  ! LHS matrix: Am c = lam c
    real(dp), intent(out) :: lam(:)   ! eigenvalues: Am c = lam c
    real(dp), intent(out) :: c(:,:)   ! eigenvectors: Am c = lam c;
                                      ! c(i,j) = ith component of jth
                                      ! vec.
    integer               :: n

    ! Lapack variables
    integer               :: lwork, liwork, info
    integer, allocatable  :: iwork(:)
    real(dp), allocatable :: work(:)
    
    preprocess: block

      ! Determine space
      n = size(Am,1)
      call assert_shape(Am, [n, n], "eigh", "Am")
      call assert_shape(c, [n, n], "eigh", "c")
      lwork = 1 + 6*n + 2*n**2
      liwork = 3 + 5*n
      allocate(work(lwork), iwork(liwork))
      c = Am

    end block preprocess

    solve: block

      ! Solve
      call dsyevd('V','L',n,c,n,lam,work,lwork,iwork,liwork,info)

    end block solve

    postprocess: block

      ! Post process
      if (info .ne. 0) then
         print *, "dsyevd returned info =", info
         if (info .lt. 0) then
            print *, "the", -info, "-th argument had an illegal value"
         else
            print *, "the algorithm failed to compute an eigenvalue while working"
            print *, "on the submatrix lying in rows and columns", 1.0_dp*info/(n+1)
            print *, "through", mod(info, n+1)
         end if
         call stop_error('eigh: dsyevd error')
      end if

    end block postprocess

  end subroutine deigh_simple

  !-------------------------------------------------------------------!
  ! Solves generalized eigen value problem for all eigenvalues and
  ! eigenvectors Am must by hermitian, Bm hermitian positive definite.
  ! Only the lower triangular part of Am and Bm is used.
  !-------------------------------------------------------------------!
  
  subroutine zeigh_generalized(Am, Bm, lam, c)

    complex(dp), intent(in)  :: Am(:,:)   ! LHS matrix: Am c = lam Bm c
    complex(dp), intent(in)  :: Bm(:,:)   ! RHS matrix: Am c = lam Bm c
    real(dp), intent(out)    :: lam(:)      ! eigenvalues: Am c = lam Bm c
    complex(dp), intent(out) :: c(:,:)   ! eigenvectors: Am c = lam Bm c; c(i,j) = ith component of jth vec.

    ! LAPACK variables
    integer                  :: info, liwork, lrwork, lwork, n
    integer, allocatable     :: iwork(:)
    real(dp), allocatable    :: rwork(:)
    complex(dp), allocatable :: Bmt(:,:), work(:)

    preprocess: block
      
      ! Determine shape
      n = size(Am,1)
      call assert_shape(Am, [n, n], "eigh", "Am")
      call assert_shape(Bm, [n, n], "eigh", "Bm")
      call assert_shape(c, [n, n], "eigh", "c")
      lwork = 2*n + n**2
      lrwork = 1 + 5*N + 2*n**2
      liwork = 3 + 5*n
      allocate(Bmt(n,n), work(lwork), rwork(lrwork), iwork(liwork))
      c = Am; Bmt = Bm  ! Bmt temporary overwritten by zhegvd

    end block preprocess
    
    solve: block
      
      ! Solve using lapack
      call zhegvd(1,'V','L',n,c,n,Bmt,n,lam,work,lwork,rwork,lrwork,iwork,liwork,info)

    end block solve

    postprocess: block

      ! Post process
      if (info .ne. 0) then
         print *, "zhegvd returned info =", info
         if (info .lt. 0) then
            print *, "the", -info, "-th argument had an illegal value"
         else if (info .le. n) then
            print *, "the algorithm failed to compute an eigenvalue while working"
            print *, "on the submatrix lying in rows and columns", 1.0_dp*info/(n+1)
            print *, "through", mod(info, n+1)
         else
            print *, "The leading minor of order ", info-n, &
                 "of B is not positive definite. The factorization of B could ", &
                 "not be completed and no eigenvalues or eigenvectors were computed."
         end if
         call stop_error('eigh: zhegvd error')
      end if
      
    end block postprocess

  end subroutine zeigh_generalized

  !-------------------------------------------------------------------!
  ! Solves eigenvalue problem for all eigenvalues and eigenvectors. Am
  ! must by symmetric Only the lower triangular part of Am is used.
  !-------------------------------------------------------------------!
  
  subroutine zeigh_simple(Am, lam, c)

    complex(dp), intent(in)  :: Am(:,:)   ! LHS matrix: Am c = lam c
    real(dp), intent(out)    :: lam(:)   ! eigenvalues: Am c = lam c
    complex(dp), intent(out) :: c(:,:)   ! eigenvectors: Am c = lam c; c(i,j) = ith component of jth vec.

    ! LAPACK variables:
    integer                  :: info, lda, liwork, lrwork, lwork, n
    integer, allocatable     :: iwork(:)
    real(dp), allocatable    :: rwork(:)
    complex(dp), allocatable :: work(:)


    preprocess: block

      ! use LAPACK's zheevd routine
      n = size(Am, 1)
      call assert_shape(Am, [n, n], "eigh", "Am")
      call assert_shape(c, [n, n], "eigh", "c")
      lda = max(1, n)
      lwork = 2*n + n**2
      lrwork = 1 + 5*n + 2*n**2
      liwork = 3 + 5*n
      allocate(work(lwork), rwork(lrwork), iwork(liwork))
      c = Am

    end block preprocess

    solve: block

      ! Solve
      call zheevd("V", "L", n, c, lda, lam, work, lwork, rwork, lrwork, &
           iwork, liwork, info)
    
    end block solve
    
    postprocess: block

      ! Post process
      if (info .ne. 0) then
         print *, "zheevd returned info =", info
         if (info .lt. 0) then
            print *, "the", -info, "-th argument had an illegal value"
         else
            print *, "the algorithm failed to compute an eigenvalue while working"
            print *, "through the submatrix lying in rows and columns through"
            print *, info/(n+1), " through ", mod(info, n+1)
         end if
         call stop_error('eigh: zheevd error')
      end if

    end block postprocess

  end subroutine zeigh_simple

  !-------------------------------------------------------------------!
  ! Return the inverse of a real matrix
  !-------------------------------------------------------------------!
  
  function dinv(Am) result(Bm)
    
    real(dp), intent(in)  :: Am(:,:)  ! matrix to be inverted
    real(dp)              :: Bm(size(Am, 1), size(Am, 2))   ! Bm = inv(Am)
    real(dp), allocatable :: Amt(:,:), work(:)  ! temporary work arrays

    ! LAPACK variables:
    integer              ::  info, lda, n, lwork, nb
    integer, allocatable :: ipiv(:)

    ! use LAPACK's dgetrf and dgetri
    n = size(Am(1, :))
    call assert_shape(Am, [n, n], "inv", "Am")
    lda = n
    nb = ilaenv(1, 'DGETRI', "UN", n, -1, -1, -1)  ! TODO: check UN param
    lwork = n*nb
    if (nb .lt. 1) nb = max(1, n)
    allocate(Amt(n,n), work(lwork), ipiv(n))
    Amt = Am

    ! Solve 
    call dgetrf(n, n, Amt, lda, ipiv, info)

    ! Post process
    if (info .ne. 0) then
       print *, "dgetrf returned info =", info
       if (info .lt. 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; The factorization"
          print *, "has been completed, but the factor U is exactly"
          print *, "singular, and division by zero will occur if it is used"
          print *, "to solve a system of equations."
       end if
       call stop_error('inv: dgetrf error')
    end if

    ! Solve
    call dgetri(n, Amt, n, ipiv, work, lwork, info)

    ! Post process
    if (info .ne. 0) then
       print *, "dgetri returned info =", info
       if (info .lt. 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; the matrix is"
          print *, "singular and its inverse could not be computed."
       end if
       call stop_error('inv: dgetri error')
    end if

    ! Result
    Bm = Amt

  end function dinv

  !-------------------------------------------------------------------!
  ! Inverse of a complex matrix
  !-------------------------------------------------------------------!
  
  function zinv(Am) result(Bm)

    ! Inverts the general complex matrix Am
    complex(dp), intent(in) :: Am(:,:)   ! Matrix to be inverted
    complex(dp)             :: Bm(size(Am, 1), size(Am, 2))   ! Bm = inv(Am)
    integer                 :: n, nb

    ! LAPACK variables
    integer                  :: lwork, info
    complex(dp), allocatable :: Amt(:,:), work(:)
    integer, allocatable     :: ipiv(:)

    ! Determine size
    n = size(Am, 1)
    call assert_shape(Am, [n, n], "inv", "Am")
    nb = ilaenv(1, 'ZGETRI', "UN", n, -1, -1, -1) ! TODO: check UN param
    if (nb .lt. 1) nb = max(1, n)
    lwork = n*nb
    allocate(Amt(n,n), ipiv(n), work(lwork))
    Amt = Am

    ! Solve
    call zgetrf(n, n, Amt, n, ipiv, info)

    ! Post process
    if (info .ne. 0) then
       print *, "zgetrf returned info =", info
       if (info .lt. 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; The factorization"
          print *, "has been completed, but the factor U is exactly"
          print *, "singular, and division by zero will occur if it is used"
          print *, "to solve a system of equations."
       end if
       call stop_error('inv: zgetrf error')
    end if

    ! Solve
    call zgetri(n, Amt, n, ipiv, work, lwork, info)
    if (info .ne. 0) then
       print *, "zgetri returned info =", info
       if (info .lt. 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; the matrix is"
          print *, "singular and its inverse could not be computed."
       end if
       call stop_error('inv: zgetri error')
    end if

    Bm = Amt

  end function zinv

  !-------------------------------------------------------------------!
  ! Solves a system of equations A x = b with one right hand side
  !-------------------------------------------------------------------!

  function dsolve(A, b) result(x)
    
    real(dp), intent(in)  :: A(:,:)  ! coefficient matrix A
    real(dp), intent(in)  :: b(:)  ! right-hand-side A x = b
    real(dp), allocatable :: x(:)

    ! LAPACK variables:
    real(dp), allocatable :: At(:,:), bt(:,:)
    integer               :: n, info, lda
    integer, allocatable  :: ipiv(:)
    
    preprocess: block

      ! Determine size
      n = size(A(1,:))
      lda = size(A(:, 1))  ! TODO: remove lda (which is = n!)
      call assert_shape(A, [n, n], "solve", "A")
      allocate(At(lda,n), bt(n,1), ipiv(n), x(n))
      At = A
      bt(:,1) = b(:)

    end block preprocess

    solve: block
      
      ! Solve the linear system
      call dgesv(n, 1, At, lda, ipiv, bt, n, info)

    end block solve

    postprocess: block

      ! Post process
      if (info .ne. 0) then
         print *, "dgesv returned info =", info
         if (info .lt. 0) then
            print *, "the", -info, "-th argument had an illegal value"
         else
            print *, "U(", info, ",", info, ") is exactly zero; The factorization"
            print *, "has been completed, but the factor U is exactly"
            print *, "singular, so the solution could not be computed."
         end if
         call stop_error('inv: dgesv error')
      end if

      ! Store the result
      x = bt(:,1)

    end block postprocess

  end function dsolve

  !-------------------------------------------------------------------!
  ! Solves a system of equations A x = b with one right hand side
  !-------------------------------------------------------------------!
  
  function zsolve(A, b) result(x)
    
    complex(dp), intent(in)  :: A(:,:)  ! coefficient matrix A
    complex(dp), intent(in)  :: b(:)    ! right-hand-side A x = b
    complex(dp), allocatable :: x(:)

    ! LAPACK variables:
    complex(dp), allocatable :: At(:,:), bt(:,:)
    integer                  :: n, info, lda
    integer, allocatable     :: ipiv(:)

    preprocess: block
      ! Determine size
      n = size(A(1,:))
      lda = size(A(:, 1))  ! TODO: remove lda here, too
      call assert_shape(A, [n, n], "solve", "A")
      allocate(At(lda,n), bt(n,1), ipiv(n), x(n))
      At = A
      bt(:,1) = b(:)
    end block preprocess
    
    solve: block
      ! Solve using LAPACK
      call zgesv(n, 1, At, lda, ipiv, bt, n, info)
    end block solve
    
    postprocess: block
      ! Post process
      if (info .ne. 0) then
         print *, "zgesv returned info =", info
         if (info .lt. 0) then
            print *, "the", -info, "-th argument had an illegal value"
         else
            print *, "U(", info, ",", info, ") is exactly zero; The factorization"
            print *, "has been completed, but the factor U is exactly"
            print *, "singular, so the solution could not be computed."
         end if
         call stop_error('inv: zgesv error')
      end if

      ! Store the result
      x = bt(:,1)
    end block postprocess

  end function zsolve
  
  !-------------------------------------------------------------------!
  ! Generates a nxn identity matrix
  !-------------------------------------------------------------------!

  pure function eye(n)

    integer, intent(in) :: n
    real(dp)    :: eye(n,n)
    integer     :: i, j

    ! generate a zero matrix
    eye = 0.0_dp
    
    ! replace the diagonals with one
    forall(i=1:n)
       eye(i,i) = 1.0_dp
    end forall

  end function eye

  !-------------------------------------------------------------------!
  ! A function that returns the determinant value of the matrix. This
  ! computes the determinant of a REAL matrix using an LU
  ! factorization.
  !-------------------------------------------------------------------!
  
  real(dp) function ddet(A) result(x)
    
    real(dp), intent(in) :: A(:, :) 
    integer              :: i

    ! LAPACK variables:
    integer               :: info, n
    integer, allocatable  :: ipiv(:)
    real(dp), allocatable :: At(:,:)

    preprocess: block
      
      ! Determine the size of the linear system
      n = size(A(1,:))
      call assert_shape(A, [n, n], "det", "A")

      !allocate(At, source=A) ! not supported by compiler yet

      ! Make a copy of the matrix as the input is modified in LU
      ! decomposition
      allocate(At(n,n), ipiv(n))
      At = A

    end block preprocess

    solve: block
      
      ! LU decomposition on the temp matrix
      call dgetrf(n, n, At, n, ipiv, info)

    end block solve

    postprocess: block
      
      ! Check the LAPACK output flag for potential issues in solution
      if (info .ne. 0) then
         print *, "dgetrf returned info =", info
         if (info .lt. 0) then
            print *, "the", -info, "-th argument had an illegal value"
         else
            print *, "U(", info, ",", info, ") is exactly zero; The factorization"
            print *, "has been completed, but the factor U is exactly"
            print *, "singular, and division by zero will occur if it is used"
            print *, "to solve a system of equations."
         end if
         call stop_error('det: dgetrf error')
      end if

      ! At now contains the LU of the factorization A = PLU as L has
      ! unit diagonal entries, the determinant can be computed from the
      ! product of U's diagonal entries. Additional sign changes
      ! stemming from the permutations P have to be taken into account
      ! as well.
      x = 1.0_dp
      do i = 1, n
         if (ipiv(i) .ne. i) then  ! additional sign change
            x = -x*At(i,i)
         else
            x = x*At(i,i)
         end if
      end do

    end block postprocess

  end function ddet
  
  !-------------------------------------------------------------------!
  ! A function that returns the determinant value of the matrix. This
  ! computes the determinant of a COMPLEX matrix using an LU
  ! factorization.
  !-------------------------------------------------------------------!
  
  complex(dp) function zdet(A) result(x)
    
    complex(dp), intent(in) :: A(:, :)
    integer                 :: i

    ! LAPACK variables:
    integer                  :: info, n
    integer, allocatable     :: ipiv(:)
    complex(dp), allocatable :: At(:,:)

    preprocess: block
      
      ! Determine the size of the linear system
      n = size(A(1,:))
      call assert_shape(A, [n, n], "det", "A")

      ! Make a copy of the matrix as the input is modified in LU
      ! decomposition
      allocate(At(n,n), ipiv(n))
      At = A

    end block preprocess

    solve: block
      
      ! LU decomposition on the temp matrix
      call zgetrf(n, n, At, n, ipiv, info)

    end block solve

    postprocess: block

      ! Check the LAPACK output flag for potential issues in solution
      if (info .ne. 0) then
         print *, "zgetrf returned info =", info
         if (info .lt. 0) then
            print *, "the", -info, "-th argument had an illegal value"
         else
            print *, "U(", info, ",", info, ") is exactly zero; The factorization"
            print *, "has been completed, but the factor U is exactly"
            print *, "singular, and division by zero will occur if it is used"
            print *, "to solve a system of equations."
         end if
         call stop_error('inv: zgetrf error')
      end if

      ! Compute the determinant (see ddet for details)
      x = 1.0_dp + 0*i_
      do i = 1, n
         if (ipiv(i) .ne. i) then  ! additional sign change
            x = -x*At(i,i)
         else
            x = x*At(i,i)
         end if
      end do

    end block postprocess

  end function zdet

  !-------------------------------------------------------------------!
  ! Compute least square solution to A x = b for real A, b
  !-------------------------------------------------------------------!

  function dlstsq(A, b) result(x)
    
    real(dp), intent(in)  :: A(:,:), b(:)
    real(dp), allocatable :: x(:)

    ! LAPACK variables:
    integer               :: info, ldb, lwork, m, n, rank
    real(dp)              :: rcond
    real(dp), allocatable :: work(:), At(:,:), Bt(:,:)
    integer, allocatable  :: jpvt(:)

    ! Handle allocations
    preprocess: block
      
      ! Determine the matrix size
      m = size(A(:,1)) ! = lda
      n = size(A(1,:))
      ldb = size(b)
      
      ! Determine the optimal workspace size
      allocate(x(n), At(m,n), Bt(ldb,1), jpvt(n), work(1))
      call dgelsy(m, n, 1, At, m, Bt, ldb, jpvt, rcond, rank, work, &
           -1, info)  ! query optimal workspace size
      lwork = int(real(work(1)))
      deallocate(work)
      
      ! Allocate with optimal size
      allocate(work(lwork))  ! allocate with ideal size
      rcond = 0.0_dp
      jpvt(:) = 0
      Bt(:,1) = b(:)  ! only one right-hand side
      At(:,:) = A(:,:)

    end block preprocess

    solve: block

      ! Call LAPACK for lease squares solution
      call dgelsy(m, n, 1, At, m, Bt, ldb, jpvt, rcond, rank, work, &
           lwork, info)

    end block solve

    postprocess: block
      
      ! Find any error code
      if (info .ne. 0) then
         print *, "dgelsy returned info = ", info
         print *, "the ", -info, "-th argument had an illegal value"
         call stop_error('lstsq: dgelsy error')
      end if
      
      ! Store the solution
      x(:) = Bt(1:n,1)

    end block postprocess

  end function dlstsq

  !-------------------------------------------------------------------!
  ! Compute least square solution to A x = b for complex A, b
  !-------------------------------------------------------------------!

  function zlstsq(A, b) result(x)

    complex(dp), intent(in)  :: A(:,:), b(:)
    complex(dp), allocatable :: x(:)

    ! LAPACK variables:
    integer                  :: info, ldb, lwork, m, n, rank
    real(dp)                 :: rcond
    complex(dp), allocatable :: At(:,:), Bt(:,:), work(:)
    real(dp), allocatable    :: rwork(:)
    integer, allocatable     :: jpvt(:)

    preprocess: block
      
      ! Determine the matrix size
      m = size(A(:,1)) ! = lda
      n = size(A(1,:))
      ldb = size(b)

      ! Query optimal workspace size
      allocate(x(n), At(m,n), Bt(ldb,1), jpvt(n), work(1), rwork(2*n))
      call zgelsy(m, n, 1, At, m, Bt, ldb, jpvt, rcond, rank, work, &
           -1, rwork, info)
      lwork = int(real(work(1)))
      deallocate(work)

      ! Allocate the optimal (ideal) workspace size
      allocate(work(lwork)) 
      rcond = 0.0_dp
      jpvt(:) = 0
      Bt(:,1) = b(:)  ! only one right-hand side
      At(:,:) = A(:,:)

    end block preprocess

    solve: block
      
      ! Call LAPACK to solve the system
      call zgelsy(m, n, 1, At, m, Bt, ldb, jpvt, rcond, rank, work, &
           lwork, rwork, info)

    end block solve
    
    ! Post process
    postprocess: block

      if (info /= 0) then
         print *, "zgelsy returned info = ", info
         print *, "the ", -info, "-th argument had an illegal value"
         call stop_error('lstsq: zgelsy error')
      end if
      
      ! Store the results
      x(:) = Bt(1:n,1)

    end block postprocess

  end function zlstsq

  !-------------------------------------------------------------------!
  ! Construct real matrix from diagonal elements
  !-------------------------------------------------------------------!
  
  pure function ddiag(x) result(A)

    real(dp), intent(in)  :: x(:)
    real(dp), allocatable :: A(:,:)
    integer               :: i, n

    n = size(x)

    allocate(A(n,n))
    A = 0.0_dp

    forall(i=1:n) A(i,i) = x(i)

  end function ddiag

  !-------------------------------------------------------------------!
  ! Construct complex matrix from diagonal elements
  !-------------------------------------------------------------------!

  pure function zdiag(x) result(A)


    complex(dp), intent(in)  :: x(:)
    complex(dp), allocatable :: A(:,:)
    integer                  :: i, n

    n = size(x)

    allocate(A(n,n))

    A = 0*i_

    forall(i=1:n) A(i,i) = x(i)

  end function zdiag

  !-------------------------------------------------------------------!
  ! Return trace along the main diagonal of a real matrix
  !-------------------------------------------------------------------!

  pure real(dp) function dtrace(A) result(t)
     
    real(dp), intent(in) :: A(:,:)
    integer              :: i

    t = 0.0_dp
    do i = 1, minval(shape(A))
       t = t + A(i,i)
    end do

  end function dtrace

  !-------------------------------------------------------------------!
  ! Return trace along the main diagonal of a real matrix
  !-------------------------------------------------------------------!

  pure complex(dp) function ztrace(A) result(t)

    complex(dp), intent(in) :: A(:,:)
    integer                 :: i

    t = 0*i_
    do i = 1, minval(shape(A))
       t = t + A(i,i)
    end do

  end function ztrace

  !-------------------------------------------------------------------!
  ! Find the singular values s_i of a real m x n matrix
  !-------------------------------------------------------------------!

  function dsvdvals(A) result(s)

    real(dp), intent(in)  :: A(:,:)
    real(dp), allocatable :: s(:)

    ! LAPACK related:
    integer               :: info, lwork, m, n
    real(dp), allocatable :: work(:), At(:,:)
    real(dp)              :: u(1,1), vt(1,1)  ! not used if only s is
                                              ! to be computed
    preprocess: block
      
      m = size(A(:,1))  ! = lda
      n = size(A(1,:))
      allocate(At(m,n), s(min(m,n)))
      At(:,:) = A(:, :)  ! A is overwritten in dgesvd

      ! Query optimal lwork and allocate workspace:   
      allocate(work(1))
      call dgesvd('N', 'N', m, n, At, m, s, u, 1, vt, 1, work, -1, info)
      lwork = int(real(work(1)))
      deallocate(work)
      allocate(work(lwork))

    end block preprocess

    solve: block
      
      ! Compute the singular value decomposition
      call dgesvd('N', 'N', m, n, At, m, s, u, 1, vt, 1, work, lwork, info)

    end block solve

    postprocess: block
      
      if (info .ne. 0) then
         print *, "dgesvd returned info = ", info
         if (info .lt. 0) then
            print *, "the ", -info, "-th argument had an illegal value"
         else
            print *, "DBDSQR did not converge, there are ", info
            print *, "superdiagonals of an intermediate bidiagonal form B"
            print *, "did not converge to zero. See the description of WORK"
            print *, "in DGESVD's man page for details."
         end if
         call stop_error('svdvals: dgesvd error')
      end if
      
    end block postprocess

  end function dsvdvals

  !-------------------------------------------------------------------!
  ! Find the singular values s_i of a complex m x n matrix
  !-------------------------------------------------------------------!
  
  function zsvdvals(A) result(s)

    complex(dp), intent(in) :: A(:,:)
    real(dp), allocatable   :: s(:)

    ! LAPACK related:
    integer                  :: info, lwork, m, n, lrwork
    complex(dp), allocatable :: work(:), At(:,:)
    real(dp), allocatable    :: rwork(:)
    complex(dp)              :: u(1,1), vt(1,1)  ! not used if only s is to be computed
    
    preprocess: block
      
      ! Determine size and allocate space
      m = size(A(:,1))  ! = lda
      n = size(A(1,:))
      lrwork = 5*min(m,n)
      allocate(At(m,n), s(min(m,n)), rwork(lrwork))
      At(:,:) = A(:,:)  ! A is overwritten in zgesvd!
      
      ! Query optimal lwork and allocate workspace
      allocate(work(1))
      call zgesvd('N', 'N', m, n, At, m, s, u, 1, vt, 1, work, -1, rwork, info)
      lwork = int(real(work(1)))
      deallocate(work)
      allocate(work(lwork))

    end block preprocess

    solve: block

      ! Compute the singular value decomposition
      call zgesvd('N', 'N', m, n, At, m, s, u, 1, vt, 1, work, lwork, rwork, info)

    end block solve

    postprocess: block
      
      if (info .ne. 0) then
         print *, "zgesvd returned info = ", info
         if (info .lt. 0) then
            print *, "the ", -info, "-th argument had an illegal value"
         else
            print *, "ZBDSQR did not converge, there are ", info
            print *, "superdiagonals of an intermediate bidiagonal form B"
            print *, "did not converge to zero. See the description of RWORK"
            print *, "in ZGESVD's man page for details."
         end if
         call stop_error('svdvals: zgesvd error')
      end if

    end block postprocess

  end function zsvdvals

  !-------------------------------------------------------------------!
  ! compute the singular value decomposition A = U sigma V^T of a real
  ! m x n matrix A. U is m x m. V^T is n x n. s has size min(m, n) -->
  ! sigma matrix is (n x m) with sigma_ii = s_i
  ! -------------------------------------------------------------------!
  
  subroutine dsvd(A, s, U, Vtransp)

    real(dp), intent(in)  :: A(:,:)
    real(dp), intent(out) :: s(:), U(:,:), Vtransp(:,:)

    ! LAPACK related:
    integer               :: info, lwork, m, n, ldu
    real(dp), allocatable :: work(:), At(:,:)

    preallocate: block

      m = size(A(:,1))  ! = lda
      n = size(A(1,:))
      ldu = m
      allocate(At(m,n))
      At(:,:) = A(:,:)  ! use a temporary as dgesvd destroys its input

      call assert_shape(U, [m, m], "svd", "U")
      call assert_shape(Vtransp, [n, n], "svd", "Vtransp")

      ! query optimal lwork and allocate workspace:
      allocate(work(1))
      call dgesvd('A', 'A', m, n, At, m, s, U, ldu, Vtransp, n, work, -1, info)
      lwork = int(real(work(1)))
      deallocate(work)
      allocate(work(lwork))

    end block preallocate

    solve: block

      call dgesvd('A', 'A', m, n, At, m, s, U, ldu, Vtransp, n, work, lwork, info)

    end block solve

    postprocess: block

      if (info /= 0) then
         print *, "dgesvd returned info = ", info
         if (info < 0) then
            print *, "the ", -info, "-th argument had an illegal value"
         else
            print *, "DBDSQR did not converge, there are ", info
            print *, "superdiagonals of an intermediate bidiagonal form B"
            print *, "did not converge to zero. See the description of WORK"
            print *, "in DGESVD's man page for details."
         end if
         call stop_error('svd: dgesvd error')
      end if

    end block postprocess

  end subroutine dsvd

  !-------------------------------------------------------------------!
  ! compute the singular value decomposition A = U sigma V^H of a
  ! complex m x m matrix A
  ! U is m x min(m, n)
  ! Vtransp is n x n
  ! sigma is m x n with with sigma_ii = s_i
  ! note that this routine returns V^H, not V!
  !-------------------------------------------------------------------!
  
  subroutine zsvd(A, s, U, Vtransp)

    complex(dp), intent(in)  :: A(:,:)
    real(dp), intent(out)    :: s(:)
    complex(dp), intent(out) :: U(:,:), Vtransp(:,:)

    ! LAPACK related:
    integer                  :: info, lwork, m, n, ldu, lrwork
    real(dp), allocatable    :: rwork(:)
    complex(dp), allocatable :: work(:), At(:,:)

    ! TODO: check shapes here and in other routines?

    preprocess: block

      m = size(A(:,1))  ! = lda
      n = size(A(1,:))
      ldu = m
      lrwork = 5*min(m,n)
      allocate(rwork(lrwork), At(m,n))
      At(:,:) = A(:,:)  ! use a temporary as zgesvd destroys its input

      call assert_shape(U, [m, m], "svd", "U")
      call assert_shape(Vtransp, [n, n], "svd", "Vtransp")

      ! query optimal lwork and allocate workspace:
      allocate(work(1))
      call zgesvd('A', 'A', m, n, At, m, s, U, ldu, Vtransp, n, work, -1,&
           rwork, info)
      lwork = int(real(work(1)))
      deallocate(work)
      allocate(work(lwork))

    end block preprocess

    solve: block

      call zgesvd('A', 'A', m, n, At, m, s, U, ldu, Vtransp, n, work, &
           lwork, rwork, info)

    end block solve

    postprocess: block

      if (info /= 0) then
         print *, "zgesvd returned info = ", info
         if (info < 0) then
            print *, "the ", -info, "-th argument had an illegal value"
         else
            print *, "ZBDSQR did not converge, there are ", info
            print *, "superdiagonals of an intermediate bidiagonal form B"
            print *, "did not converge to zero. See the description of WORK"
            print *, "in DGESVD's man page for details."
         end if
         call stop_error('svd: zgesvd error')
      end if

    end block postprocess

  end subroutine zsvd
  
  !-------------------------------------------------------------------!
  ! Make sure a given real matrix has a given shape
  !-------------------------------------------------------------------!

  subroutine dassert_shape(A, shap, routine, matname)

    real(dp), intent(in) :: A(:,:)
    integer, intent(in)  :: shap(:)
    character(len=*)     :: routine, matname

    if (any(shape(A) .ne. shap)) then
       print *, "In routine " // routine // " matrix " // matname // " has illegal shape ", shape(A)
       print *, "Shape should be ", shap
       call stop_error("Aborting due to illegal matrix operation")
    end if

  end subroutine dassert_shape
 
  !-------------------------------------------------------------------!
  ! Make sure a given real matrix has a given shape
  !-------------------------------------------------------------------!
  
  subroutine zassert_shape(A, shap, routine, matname)
    
    complex(dp), intent(in) :: A(:,:)
    integer, intent(in)     :: shap(:)
    character(len=*)        :: routine, matname
    
    if (any(shape(A) .ne. shap)) then
       print *, "In routine " // routine // " matrix " // matname // " has illegal shape ", shape(A)
       print *, "Shape should be ", shap
       call stop_error("Aborting due to illegal matrix operation")
    end if

  end subroutine zassert_shape
  
  !-------------------------------------------------------------------!
  ! Stop with an error message
  !-------------------------------------------------------------------!

  subroutine stop_error( message )

    character(len=*) :: message

    print *, message
    stop    

  end subroutine stop_error

  !-------------------------------------------------------------------!
  ! Assert if the condition is not met (i.e. condition .eq. false)
  !-------------------------------------------------------------------!

  subroutine assert(condition)

    logical, intent(in) :: condition

    if (.not. condition) call stop_error("Assert failed.")

  end subroutine assert

end module linear_algebra
