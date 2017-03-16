!=====================================================================!
! Module that wraps important LAPACK calls. The implementation is not
! provided and is expected to be in the form of compiled library or
! overridden procedures by the user.
!=====================================================================!

module lapack
  
  use iso_fortran_env, only : dp => REAL64

  implicit none

  ! This is the precision that LAPACK "d" routines were compiled with
  ! (typically double precision, unless a special compiler option was
  ! used while compiling LAPACK). This "dp" is only used in lapack.f90
  ! The "d" routines data type is defined as "double precision", so we
  ! make "dp" the same kind as 0.d0 ("double precision"), so as long
  ! as LAPACK and this file were compiled with the same compiler
  ! options, it will be consistent. (If for example all double
  ! precision is promoted to quadruple precision, it will be promoted
  ! both in LAPACK and here.)

  interface

     SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       import :: dp
       INTEGER            INFO, LDA, LDB, N, NRHS
       INTEGER            IPIV( * )
       REAL(dp)           A( LDA, * ), B( LDB, * )
     END SUBROUTINE DGESV

     SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, &
          EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, &
          IWORK, INFO )
       import :: dp
       CHARACTER          EQUED, FACT, TRANS
       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
       REAL(dp)           RCOND
       INTEGER            IPIV( * ), IWORK( * )
       REAL(dp)           A( LDA, * ), AF( LDAF, * ), B( LDB, * ), BERR( * ), &
            C( * ), FERR( * ), R( * ), WORK( * ), X( LDX, * )
     END SUBROUTINE DGESVX

     SUBROUTINE ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
       import :: dp
       INTEGER            INFO, LDA, LDB, N, NRHS
       INTEGER            IPIV( * )
       COMPLEX(dp)        A( LDA, * ), B( LDB, * )
     END SUBROUTINE ZGESV

     SUBROUTINE ZGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, &
          EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, &
          WORK, RWORK, INFO )
       import :: dp
       CHARACTER          EQUED, FACT, TRANS
       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
       REAL(dp)           RCOND
       INTEGER            IPIV( * )
       REAL(dp)           BERR( * ), C( * ), FERR( * ), R( * ), RWORK( * )
       COMPLEX(dp)        A( LDA, * ), AF( LDAF, * ), B( LDB, * ), WORK( * ), &
            X( LDX, * )
     END SUBROUTINE ZGESVX

     SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
       import :: dp
       INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
       INTEGER            IPIV( * )
       REAL(dp)           AB( LDAB, * ), B( LDB, * )
     END SUBROUTINE DGBSV

     SUBROUTINE DSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
       import :: dp
       CHARACTER          UPLO
       INTEGER            INFO, LDA, LDB, LWORK, N, NRHS
       INTEGER            IPIV( * )
       REAL(dp)           A( LDA, * ), B( LDB, * ), WORK( * )
     END SUBROUTINE DSYSV

     SUBROUTINE DSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, &
          LDB, X, LDX, RCOND, FERR, BERR, WORK, LWORK, &
          IWORK, INFO )
       import :: dp
       CHARACTER          FACT, UPLO
       INTEGER            INFO, LDA, LDAF, LDB, LDX, LWORK, N, NRHS
       REAL(dp)           RCOND
       INTEGER            IPIV( * ), IWORK( * )
       REAL(dp)           A( LDA, * ), AF( LDAF, * ), B( LDB, * ), &
            BERR( * ), FERR( * ), WORK( * ), X( LDX, * )
     END SUBROUTINE DSYSVX

     SUBROUTINE DSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, &
          LIWORK, INFO )
       import :: dp
       CHARACTER          JOBZ, UPLO
       INTEGER            INFO, LDA, LIWORK, LWORK, N
       INTEGER            IWORK( * )
       REAL(dp)           A( LDA, * ), W( * ), WORK( * )
     END SUBROUTINE DSYEVD

     SUBROUTINE DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB, &
          VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, &
          LWORK, IWORK, IFAIL, INFO )
       import :: dp
       CHARACTER          JOBZ, RANGE, UPLO
       INTEGER            IL, INFO, ITYPE, IU, LDA, LDB, LDZ, LWORK, M, N
       REAL(dp)           ABSTOL, VL, VU
       INTEGER            IFAIL( * ), IWORK( * )
       REAL(dp)           A( LDA, * ), B( LDB, * ), W( * ), WORK( * ), &
            Z( LDZ, * )
     END SUBROUTINE DSYGVX

     SUBROUTINE DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, &
          BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
       import :: dp
       CHARACTER          JOBVL, JOBVR
       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
       REAL(dp)           A( LDA, * ), ALPHAI( * ), ALPHAR( * ), &
            B( LDB, * ), BETA( * ), VL( LDVL, * ), &
            VR( LDVR, * ), WORK( * )
     END SUBROUTINE DGGEV

     SUBROUTINE DGGEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, B, LDB, &
          ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, ILO, IHI, &
          LSCALE, RSCALE, ABNRM, BBNRM, RCONDE, RCONDV, WORK, &
          LWORK, IWORK, BWORK, INFO )
       import :: dp
       CHARACTER          BALANC, JOBVL, JOBVR, SENSE
       INTEGER            IHI, ILO, INFO, LDA, LDB, LDVL, LDVR, LWORK, N
       REAL(dp)           ABNRM, BBNRM
       LOGICAL            BWORK( * )
       INTEGER            IWORK( * )
       REAL(dp)           A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), &
            BETA( * ), LSCALE( * ), RCONDE( * ), RCONDV( * ), &
            RSCALE( * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * )
     END SUBROUTINE DGGEVX

     SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, &
          LDVR, WORK, LWORK, INFO )
       import :: dp
       CHARACTER          JOBVL, JOBVR
       INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
       REAL(dp)           A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), WI( * ), &
            WORK( * ), WR( * )
     END SUBROUTINE DGEEV

     SUBROUTINE DGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR, WI, &
          VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, &
          RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )
       import :: dp
       CHARACTER          BALANC, JOBVL, JOBVR, SENSE
       INTEGER            IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N
       REAL(dp)           ABNRM
       INTEGER            IWORK( * )
       REAL(dp)           A( LDA, * ), RCONDE( * ), RCONDV( * ), &
            SCALE( * ), VL( LDVL, * ), VR( LDVR, * ), &
            WI( * ), WORK( * ), WR( * )
     END SUBROUTINE DGEEVX

     SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, &
          WORK, LWORK, RWORK, INFO )
       import :: dp
       CHARACTER          JOBVL, JOBVR
       INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
       REAL(dp)           RWORK( * )
       COMPLEX(dp)        A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), W( * ), &
            WORK( * )
     END SUBROUTINE ZGEEV

     SUBROUTINE ZGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, W, VL, &
          LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, &
          RCONDV, WORK, LWORK, RWORK, INFO )
       import :: dp
       CHARACTER          BALANC, JOBVL, JOBVR, SENSE
       INTEGER            IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N
       REAL(dp)           ABNRM
       REAL(dp)           RCONDE( * ), RCONDV( * ), RWORK( * ), SCALE( * )
       COMPLEX(dp)        A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), W( * ), &
            WORK( * )
     END SUBROUTINE ZGEEVX

     SUBROUTINE DSYGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, &
          LWORK, IWORK, LIWORK, INFO )
       import :: dp
       CHARACTER          JOBZ, UPLO
       INTEGER            INFO, ITYPE, LDA, LDB, LIWORK, LWORK, N
       INTEGER            IWORK( * )
       REAL(dp)           A( LDA, * ), B( LDB, * ), W( * ), WORK( * )
     END SUBROUTINE DSYGVD

     REAL(dp) FUNCTION DLAMCH( CMACH )
       import :: dp
       CHARACTER          CMACH
     END FUNCTION DLAMCH

     INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
       CHARACTER*( * )    NAME, OPTS
       INTEGER            ISPEC, N1, N2, N3, N4
     END FUNCTION ILAENV

     SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
       import :: dp
       INTEGER            INFO, LDA, M, N
       INTEGER            IPIV( * )
       COMPLEX(dp)        A( LDA, * )
     END SUBROUTINE ZGETRF

     SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
       import :: dp
       INTEGER            INFO, LDA, LWORK, N
       INTEGER            IPIV( * )
       COMPLEX(dp)        A( LDA, * ), WORK( * )
     END SUBROUTINE ZGETRI

     SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
       import :: dp
       INTEGER            INFO, LDA, M, N
       INTEGER            IPIV( * )
       REAL(dp)           A( LDA, * )
     END SUBROUTINE DGETRF

     SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
       import :: dp
       INTEGER            INFO, LDA, LWORK, N
       INTEGER            IPIV( * )
       REAL(dp)           A( LDA, * ), WORK( * )
     END SUBROUTINE DGETRI

     SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )
       import :: dp
       CHARACTER          JOBZ, UPLO
       INTEGER            INFO, LDA, LWORK, N
       REAL(dp)           RWORK( * ), W( * )
       COMPLEX(dp)        A( LDA, * ), WORK( * )
     END SUBROUTINE ZHEEV

     SUBROUTINE ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, &
          LRWORK, IWORK, LIWORK, INFO )
       import :: dp
       CHARACTER          JOBZ, UPLO
       INTEGER            INFO, LDA, LIWORK, LRWORK, LWORK, N
       INTEGER            IWORK( * )
       REAL(dp)           RWORK( * ), W( * )
       COMPLEX(dp)        A( LDA, * ), WORK( * )
     END SUBROUTINE ZHEEVD

     SUBROUTINE ZHEGVD( ITYPE,  JOBZ,  UPLO,  N,  A,  LDA,  B, LDB, W, &
          WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, &
          INFO )
       import :: dp
       CHARACTER          JOBZ, UPLO
       INTEGER            INFO, ITYPE, LDA, LDB, LIWORK, LRWORK, LWORK, N
       INTEGER            IWORK( * )
       REAL(dp)           RWORK( * ), W( * )
       COMPLEX(dp)        A( LDA, * ), B( LDB, * ), WORK( * )
     END SUBROUTINE ZHEGVD

     SUBROUTINE DGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, &
          WORK, LWORK, INFO )
       import :: dp
       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
       REAL(dp)           RCOND
       INTEGER            JPVT( * )
       REAL(dp)           A( LDA, * ), B( LDB, * ), WORK( * )
     END SUBROUTINE DGELSY

     SUBROUTINE ZGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, &
          WORK, LWORK, RWORK, INFO )
       import :: dp
       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
       REAL(dp)           RCOND
       INTEGER            JPVT( * )
       REAL(dp)           RWORK( * )
       COMPLEX(dp)        A( LDA, * ), B( LDB, * ), WORK( * )
     END SUBROUTINE ZGELSY

     SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, &
          LDVT, WORK, LWORK, INFO )
       import :: dp
       CHARACTER          JOBU, JOBVT
       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
       REAL(dp)           A( LDA, * ), S( * ),  U( LDU,  * ), VT( LDVT, * ), &
            WORK( * )
     END SUBROUTINE DGESVD

     SUBROUTINE ZGESVD( JOBU, JOBVT,  M,  N,  A,  LDA, S, U, LDU, VT, LDVT, &
          WORK, LWORK, RWORK, INFO )
       import :: dp
       CHARACTER          JOBU, JOBVT
       INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
       REAL(dp)           RWORK( * ), S( * )
       COMPLEX(dp)        A( LDA, * ), U( LDU, * ), VT( LDVT, * ), WORK( * )
     END SUBROUTINE ZGESVD

  end interface

contains

end module lapack
