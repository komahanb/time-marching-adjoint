!=====================================================================!
! Interface module for BLAS routines. This code is based on
! Ondrej Certik (https://github.com/certik) LANL, NM. 
!=====================================================================!

module blas

  use constants, only : WP

  implicit none

  ! This is the precision that BLAS "d" routines were compiled with
  ! (typically double precision, unless a special compiler option was
  ! used while compiling BLAS). This "dp" is only used in lapack.f90
  ! The "d" routines data type is defined as "double precision", so we
  ! make "dp" the same kind as 0._WP ("double precision"), so as long
  ! as BLAS and this file were compiled with the same compiler
  ! options, it will be consistent. (If for example all double
  ! precision is promoted to quadruple precision, it will be promoted
  ! both in BLAS and here.)

  interface

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

end module blas
