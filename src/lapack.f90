!=====================================================================!
! Interface module for LAPACK routines. This code is based on Ondrej
! Certik (https://github.com/certik) LANL, NM.
!
! See http://www.netlib.org/lapack/double/ and complex for more
! details.
! =====================================================================!

module lapack

  use constants, only : WP

  implicit none

  interface

     !================================================================!
     ! LINEAR SYSTEM SOLUTION ROUTINES FOLLOW
     !================================================================!

     !----------------------------------------------------------------!
     ! Solves a general system of linear equations AX=B
     !----------------------------------------------------------------!
     
     SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       import :: WP
       INTEGER            INFO, LDA, LDB, N, NRHS
       INTEGER            IPIV( * )
       REAL(WP)           A( LDA, * ), B( LDB, * )
     END SUBROUTINE DGESV

     !----------------------------------------------------------------!
     ! Solves a general system of linear equations AX=B, A**T X=B or
     ! A**H X=B, and provides an estimate of the condition number and
     ! error bounds on the solution.
     ! ----------------------------------------------------------------!
     
     SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, &
          EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, &
          WORK, IWORK, INFO )
       import :: WP
       CHARACTER          EQUED, FACT, TRANS
       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
       REAL(WP)           RCOND
       INTEGER            IPIV( * ), IWORK( * )
       REAL(WP)           A( LDA, * ), AF( LDAF, * ), B( LDB, * ), BERR( * ), &
            C( * ), FERR( * ), R( * ), WORK( * ), X( LDX, * )
     END SUBROUTINE DGESVX
     
     !----------------------------------------------------------------!
     ! Solves a general system of linear equations AX=B
     !----------------------------------------------------------------!
     
     SUBROUTINE ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       import :: WP
       INTEGER            INFO, LDA, LDB, N, NRHS
       INTEGER            IPIV( * )
       COMPLEX(WP)        A( LDA, * ), B( LDB, * )
     END SUBROUTINE ZGESV
     
     !----------------------------------------------------------------!
     ! Solves a general system of linear equations AX=B, A**T X=B or
     ! A**H X=B, and provides an estimate of the condition number and
     ! error bounds on the solution.
     ! ----------------------------------------------------------------!

     SUBROUTINE ZGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, &
          EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, &
          WORK, RWORK, INFO )
       import :: WP
       CHARACTER          EQUED, FACT, TRANS
       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
       REAL(WP)           RCOND
       INTEGER            IPIV( * )
       REAL(WP)           BERR( * ), C( * ), FERR( * ), R( * ), RWORK( * )
       COMPLEX(WP)        A( LDA, * ), AF( LDAF, * ), B( LDB, * ), WORK( * ), &
            X( LDX, * )
     END SUBROUTINE ZGESVX
     
     !-----------------------------------------------------------------!
     ! Solves a general banded system of linear equations AX=B.
     !-----------------------------------------------------------------!
     
     SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
       import :: WP
       INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
       INTEGER            IPIV( * )
       REAL(WP)           AB( LDAB, * ), B( LDB, * )
     END SUBROUTINE DGBSV
     
     !------------------------------------------------------------------!
     ! Solves a real symmetric indefinite system of linear equations
     ! AX=B.
     !------------------------------------------------------------------!
     
     SUBROUTINE DSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
       import :: WP
       CHARACTER          UPLO
       INTEGER            INFO, LDA, LDB, LWORK, N, NRHS
       INTEGER            IPIV( * )
       REAL(WP)           A( LDA, * ), B( LDB, * ), WORK( * )
     END SUBROUTINE DSYSV

     !------------------------------------------------------------------!
     ! Solves a real symmetric indefinite system of linear equations
     ! AX=B, and provides an estimate of the condition number and
     ! error bounds on the solution.
     !-------------------------------------------------------------------!
     
     SUBROUTINE DSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, &
          LDB, X, LDX, RCOND, FERR, BERR, WORK, LWORK, &
          IWORK, INFO )
       import :: WP
       CHARACTER          FACT, UPLO
       INTEGER            INFO, LDA, LDAF, LDB, LDX, LWORK, N, NRHS
       REAL(WP)           RCOND
       INTEGER            IPIV( * ), IWORK( * )
       REAL(WP)           A( LDA, * ), AF( LDAF, * ), B( LDB, * ), &
            BERR( * ), FERR( * ), WORK( * ), X( LDX, * )
     END SUBROUTINE DSYSVX

     !==================================================================!
     ! EIGEN VECTOR/ VALUE ROUTINES FOLLOW
     !==================================================================!
     
     !------------------------------------------------------------------!
     ! Computes all eigenvalues, and optionally, eigenvectors of a
     ! real symmetric matrix.  If eigenvectors are desired, it uses a
     ! divide and conquer algorithm.
     ! ------------------------------------------------------------------!
     
     SUBROUTINE DSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, &
          LIWORK, INFO )
       import :: WP
       CHARACTER          JOBZ, UPLO
       INTEGER            INFO, LDA, LIWORK, LWORK, N
       INTEGER            IWORK( * )
       REAL(WP)           A( LDA, * ), W( * ), WORK( * )
     END SUBROUTINE DSYEVD

     !------------------------------------------------------------------!
     ! Computes selected eigenvalues, and optionally, the eigenvectors
     ! of a generalized symmetric-definite generalized eigenproblem,
     ! Ax= lambda Bx, ABx= lambda x, or BAx= lambda x.
     ! ------------------------------------------------------------------!
     
     SUBROUTINE DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB, &
          VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, &
          LWORK, IWORK, IFAIL, INFO )
       import :: WP
       CHARACTER          JOBZ, RANGE, UPLO
       INTEGER            IL, INFO, ITYPE, IU, LDA, LDB, LDZ, LWORK, M, N
       REAL(WP)           ABSTOL, VL, VU
       INTEGER            IFAIL( * ), IWORK( * )
       REAL(WP)           A( LDA, * ), B( LDB, * ), W( * ), WORK( * ), &
            Z( LDZ, * )
     END SUBROUTINE DSYGVX

     !-------------------------------------------------------------------!
     ! Computes the generalized eigenvalues, and optionally, the left
     ! and/or right generalized eigenvectors.
     !-------------------------------------------------------------------!
     
     SUBROUTINE DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, &
          BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
       import :: WP
       CHARACTER          JOBVL, JOBVR
       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
       REAL(WP)           A( LDA, * ), ALPHAI( * ), ALPHAR( * ), &
            B( LDB, * ), BETA( * ), VL( LDVL, * ), &
            VR( LDVR, * ), WORK( * )
     END SUBROUTINE DGGEV

     !-------------------------------------------------------------------!
     ! Computes the generalized eigenvalues, and optionally, the left
     ! and/or right generalized eigenvectors.
     ! -------------------------------------------------------------------!
     
     SUBROUTINE DGGEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, B, LDB, &
          ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, ILO, IHI, &
          LSCALE, RSCALE, ABNRM, BBNRM, RCONDE, RCONDV, WORK, &
          LWORK, IWORK, BWORK, INFO )
       import :: WP
       CHARACTER          BALANC, JOBVL, JOBVR, SENSE
       INTEGER            IHI, ILO, INFO, LDA, LDB, LDVL, LDVR, LWORK, N
       REAL(WP)           ABNRM, BBNRM
       LOGICAL            BWORK( * )
       INTEGER            IWORK( * )
       REAL(WP)           A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), &
            BETA( * ), LSCALE( * ), RCONDE( * ), RCONDV( * ), &
            RSCALE( * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * )
     END SUBROUTINE DGGEVX
     
     !-------------------------------------------------------------------!
     ! Computes the eigenvalues and left and right eigenvectors of a
     ! general matrix, with preliminary balancing of the matrix, and
     ! computes reciprocal condition numbers for the eigenvalues and
     ! right eigenvectors.
     !-------------------------------------------------------------------!
     
     SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, &
          LDVR, WORK, LWORK, INFO )
       import :: WP
       CHARACTER          JOBVL, JOBVR
       INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
       REAL(WP)           A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), WI( * ), &
            WORK( * ), WR( * )
     END SUBROUTINE DGEEV

     !-------------------------------------------------------------------!
     ! Computes the generalized eigenvalues, and optionally, the left
     ! and/or right generalized eigenvectors.
     !-------------------------------------------------------------------!
          
     SUBROUTINE DGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR, WI, &
          VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, &
          RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )
       import :: WP
       CHARACTER          BALANC, JOBVL, JOBVR, SENSE
       INTEGER            IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N
       REAL(WP)           ABNRM
       INTEGER            IWORK( * )
       REAL(WP)           A( LDA, * ), RCONDE( * ), RCONDV( * ), &
            SCALE( * ), VL( LDVL, * ), VR( LDVR, * ), &
            WI( * ), WORK( * ), WR( * )
     END SUBROUTINE DGEEVX
     

     SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, &
          WORK, LWORK, RWORK, INFO )
       import :: WP
       CHARACTER          JOBVL, JOBVR
       INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
       REAL(WP)           RWORK( * )
       COMPLEX(WP)        A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), W( * ), &
            WORK( * )
     END SUBROUTINE ZGEEV

     SUBROUTINE ZGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, W, VL, &
          LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, &
          RCONDV, WORK, LWORK, RWORK, INFO )
       import :: WP
       CHARACTER          BALANC, JOBVL, JOBVR, SENSE
       INTEGER            IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N
       REAL(WP)           ABNRM
       REAL(WP)           RCONDE( * ), RCONDV( * ), RWORK( * ), SCALE( * )
       COMPLEX(WP)        A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), W( * ), &
            WORK( * )
     END SUBROUTINE ZGEEVX

     !----------------------------------------------------------------------!
     !  Computes all eigenvalues and the eigenvectors of a generalized
     !  symmetric-definite generalized eigenproblem, Ax= lambda Bx,
     !  ABx= lambda x, or BAx= lambda x.  If eigenvectors are desired,
     !  it uses a divide and conquer algorithm.
     !----------------------------------------------------------------------!

     SUBROUTINE DSYGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, &
          LWORK, IWORK, LIWORK, INFO )
       import :: WP
       CHARACTER          JOBZ, UPLO
       INTEGER            INFO, ITYPE, LDA, LDB, LIWORK, LWORK, N
       INTEGER            IWORK( * )
       REAL(WP)           A( LDA, * ), B( LDB, * ), W( * ), WORK( * )
     END SUBROUTINE DSYGVD

     !----------------------------------------------------------------------!
     ! Determines double precision machine parameters.
     !----------------------------------------------------------------------!
     
     REAL(WP) FUNCTION DLAMCH( CMACH )
       import :: WP
       CHARACTER          CMACH
     END FUNCTION DLAMCH

     INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
       CHARACTER*( * )    NAME, OPTS
       INTEGER            ISPEC, N1, N2, N3, N4
     END FUNCTION ILAENV

     SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
       import :: WP
       INTEGER            INFO, LDA, M, N
       INTEGER            IPIV( * )
       COMPLEX(WP)        A( LDA, * )
     END SUBROUTINE ZGETRF

     SUBROUTINE ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       import :: WP
       CHARACTER          TRANS
       INTEGER            INFO, LDA, LDB, N, NRHS
       INTEGER            IPIV( * )
       COMPLEX(WP)         A( LDA, * ), B( LDB, * )
     END SUBROUTINE ZGETRS

     SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
       import :: WP
       INTEGER            INFO, LDA, LWORK, N
       INTEGER            IPIV( * )
       COMPLEX(WP)        A( LDA, * ), WORK( * )
     END SUBROUTINE ZGETRI

     !----------------------------------------------------------------------!
     ! Computes an LU factorization of a general matrix, using partial
     ! pivoting with row interchanges.
     ! ----------------------------------------------------------------------!
     
     SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
       import :: WP
       INTEGER            INFO, LDA, M, N
       INTEGER            IPIV( * )
       REAL(WP)           A( LDA, * )
     END SUBROUTINE DGETRF

     !----------------------------------------------------------------------!
     ! Computes the inverse of a general matrix, using the LU
     ! factorization computed by DGETRF.
     ! ----------------------------------------------------------------------!
     
     SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
       import :: WP
       INTEGER            INFO, LDA, LWORK, N
       INTEGER            IPIV( * )
       REAL(WP)           A( LDA, * ), WORK( * )
     END SUBROUTINE DGETRI

     SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )
       import :: WP
       CHARACTER          JOBZ, UPLO
       INTEGER            INFO, LDA, LWORK, N
       REAL(WP)           RWORK( * ), W( * )
       COMPLEX(WP)        A( LDA, * ), WORK( * )
     END SUBROUTINE ZHEEV

     SUBROUTINE ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, &
          LRWORK, IWORK, LIWORK, INFO )
       import :: WP
       CHARACTER          JOBZ, UPLO
       INTEGER            INFO, LDA, LIWORK, LRWORK, LWORK, N
       INTEGER            IWORK( * )
       REAL(WP)           RWORK( * ), W( * )
       COMPLEX(WP)        A( LDA, * ), WORK( * )
     END SUBROUTINE ZHEEVD

     SUBROUTINE ZHEGVD( ITYPE,  JOBZ,  UPLO,  N,  A,  LDA,  B, LDB, W, &
          WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, &
          INFO )
       import :: WP
       CHARACTER          JOBZ, UPLO
       INTEGER            INFO, ITYPE, LDA, LDB, LIWORK, LRWORK, LWORK, N
       INTEGER            IWORK( * )
       REAL(WP)           RWORK( * ), W( * )
       COMPLEX(WP)        A( LDA, * ), B( LDB, * ), WORK( * )
     END SUBROUTINE ZHEGVD
     
     !----------------------------------------------------------------------!
     ! Computes the minimum norm least squares solution to an over- or
     ! under-determined system of linear equations A X=B, using a
     ! complete orthogonal factorization of A.
     ! ----------------------------------------------------------------------!
     
     SUBROUTINE DGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, &
          WORK, LWORK, INFO )
       import :: WP
       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
       REAL(WP)           RCOND
       INTEGER            JPVT( * )
       REAL(WP)           A( LDA, * ), B( LDB, * ), WORK( * )
     END SUBROUTINE DGELSY

     SUBROUTINE ZGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, &
          WORK, LWORK, RWORK, INFO )
       import :: WP
       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
       REAL(WP)           RCOND
       INTEGER            JPVT( * )
       REAL(WP)           RWORK( * )
       COMPLEX(WP)        A( LDA, * ), B( LDB, * ), WORK( * )
     END SUBROUTINE ZGELSY

     !----------------------------------------------------------------------!
     ! Computes the singular value decomposition (SVD) of a general
     ! rectangular matrix.
     !----------------------------------------------------------------------!
     
     SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, &
          LDVT, WORK, LWORK, INFO )
       import :: WP
       CHARACTER          JOBU, JOBVT
       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
       REAL(WP)           A( LDA, * ), S( * ),  U( LDU,  * ), VT( LDVT, * ), &
            WORK( * )
     END SUBROUTINE DGESVD

     SUBROUTINE ZGESVD( JOBU, JOBVT,  M,  N,  A,  LDA, S, U, LDU, VT, LDVT, &
          WORK, LWORK, RWORK, INFO )
       import :: WP
       CHARACTER          JOBU, JOBVT
       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
       REAL(WP)           RWORK( * ), S( * )
       COMPLEX(WP)        A( LDA, * ), U( LDU, * ), VT( LDVT, * ), WORK( * )
     END SUBROUTINE ZGESVD

     SUBROUTINE XERBLA( SRNAME, INFO )
       CHARACTER*(*)      SRNAME
       INTEGER            INFO
     END SUBROUTINE XERBLA

     !================================================================!
     ! BLAS
     !================================================================!

     SUBROUTINE ZCOPY(N,ZX,INCX,ZY,INCY)
       import :: WP
       INTEGER INCX,INCY,N
       COMPLEX(WP) ZX(*),ZY(*)
     END SUBROUTINE ZCOPY

     SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
       import :: WP
       DOUBLE PRECISION ALPHA,BETA
       INTEGER K,LDA,LDB,LDC,M,N
       CHARACTER TRANSA,TRANSB
       REAL(WP) A(LDA,*),B(LDB,*),C(LDC,*)
     END SUBROUTINE DGEMM

     SUBROUTINE DSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
       import :: WP
       REAL(WP) ALPHA,BETA
       INTEGER LDA,LDB,LDC,M,N
       CHARACTER SIDE,UPLO
       REAL(WP) A(LDA,*),B(LDB,*),C(LDC,*)
     END SUBROUTINE DSYMM

  end interface

contains

end module lapack
