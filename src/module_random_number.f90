!=====================================================================!
! Interface for generating uniform random numbers
! 
! Author : Komahan Boopathy
!=====================================================================!

module interface_random

  use iso_fortran_env, only : dp => REAL64

  implicit none
  
  type, abstract :: random

   contains

     ! Deferred procedure
     procedure(get_interface), deferred :: get

     ! scale between different bounds
     procedure :: scale

  end type random
  
  interface

     !================================================================!
     ! Interface for getting random numbers of different kind
     !================================================================!

     impure elemental subroutine get_interface(this, x, xlow, xup)

       import random
       import dp

       class(random) , intent(in)  :: this 
       real(dp)      , intent(out) :: x
       real(dp)      , intent(in), optional :: xlow
       real(dp)      , intent(in), optional :: xup    

     end subroutine get_interface

  end interface

contains

  !===================================================================!
  ! Offset the point between supplied lower and upper bounds
  !===================================================================!
  
  pure elemental subroutine scale(this, x, xlow, xup)

    class(random) , intent(in)    :: this 
    real(dp)      , intent(inout) :: x
    real(dp)      , intent(in)    :: xlow
    real(dp)      , intent(in)    :: xup

    x = xlow + (xup - xlow)*x

  end subroutine scale

end module interface_random

!=====================================================================!
! Class for returning uniform random numbers
! 
! Author : Komahan Boopathy
!=====================================================================!

module class_random_uniform
  
  use iso_fortran_env , only : dp => REAL64  
  use interface_random, only : random

  implicit none

  ! Type for generating random numbers following uniform distribution
  type, extends(random) :: uniform 
   contains
     procedure :: get
  end type uniform
  
contains
  
  !===================================================================!
  ! Specific implementation for returning uniform random numbers
  !===================================================================!
  
  impure elemental subroutine get(this, x, xlow, xup)

    class(uniform), intent(in) :: this 
    real(dp)      , intent(out) :: x
    real(dp)      , intent(in), optional :: xlow
    real(dp)      , intent(in), optional :: xup    

    ! Generate a uniform random number using intrinsic
    call random_number(x)

    ! Scale the numbers between lower and upper bound if sought. Else
    ! the number is between [0,1]
    if (present(xlow) .and. present(xup)) then
       call this % scale(x, xlow, xup)
    end if

  end subroutine get
  
end module class_random_uniform

!=====================================================================!
! Class for returning normal random numbers
! 
! Author : Komahan Boopathy
!=====================================================================!

module class_random_normal
  
  use iso_fortran_env , only : dp => REAL64  
  use interface_random, only : random

  implicit none

  ! Type for generating random numbers following normal distribution
  type, extends(random) :: normal 
     real(dp) :: mu
     real(dp) :: sigma
   contains
     procedure :: get
     procedure :: inverse_normal
  end type normal
  
  interface normal
     module procedure create_normal
  end interface normal

contains
  
  !===================================================================!
  ! Constructor for normal distribution
  !===================================================================!
  
  pure type(normal) function create_normal(mu, sigma) result (this)

    real(dp), intent(in) :: mu, sigma

    this % mu = mu
    this % sigma = sigma

  end function create_normal

  !===================================================================!
  ! Specific implementation for returning normal random numbers
  !===================================================================!
  
  impure elemental subroutine get(this, x, xlow, xup)

    class(normal) , intent(in)  :: this 
    real(dp)      , intent(out) :: x
    real(dp)      , intent(in), optional :: xlow
    real(dp)      , intent(in), optional :: xup    

    ! Generate a normal random number using intrinsic
    call random_number(x)

    ! Distribute it normally (standard distribution)
    call this % inverse_normal(x)
    
    ! Scale this for other distributions
    x = this % mu + this % sigma * x

    ! Scale the numbers between lower and upper bound if sought. Else
    ! the number is between [0,1]
!!$    if (present(xlow) .and. present(xup)) then
!!$       call this % scale(x, xlow, xup)
!!$    end if

  end subroutine get
  
  pure elemental subroutine inverse_normal(this, p)

    class(normal) , intent(in)    :: this         
    real(dp)      , intent(inout) :: p

    real(dp) :: p_low, p_high
    real(dp) :: a1, a2, a3, a4, a5, a6
    real(dp) :: b1, b2, b3, b4, b5
    real(dp) :: c1, c2, c3, c4, c5, c6
    real(dp) :: d1, d2, d3, d4
    real(dp) :: z, q, r

    a1 = -39.6968302866538
    a2 = 220.946098424521
    a3 = -275.928510446969
    a4 = 138.357751867269
    a5 = -30.6647980661472
    a6 = 2.50662827745924

    b1 = -54.4760987982241
    b2 = 161.585836858041
    b3 = -155.698979859887
    b4 = 66.8013118877197
    b5 = -13.2806815528857

    c1 = -0.00778489400243029
    c2 = -0.322396458041136
    c3 = -2.40075827716184
    c4 = -2.54973253934373
    c5 = 4.37466414146497
    c6 = 2.93816398269878

    d1 = 0.00778469570904146
    d2 = 0.32246712907004
    d3 = 2.445134137143
    d4 = 3.75440866190742

    p_low  = 0.02425
    p_high = 1.0_dp - p_low

    if (p .lt. p_low) goto 201
    if (p .ge. p_low) goto 301

201 q = sqrt(-2.0_dp*log(p))
    z = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)
    goto 204

301 if ((p .ge. p_low) .and. (p .le. p_high)) goto 202
    if (p .gt. p_high) goto 302

202 q = p - 0.5_dp
    r = q*q
    z = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
    goto 204

302 if ((p.gt.p_high) .and. (p.lt.1)) goto 203
    
203 q = sqrt(-2*dlog(1-p))
    z = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1)

204 p = z

  end subroutine inverse_normal
  
end module class_random_normal

subroutine test

  use iso_fortran_env     , only : dp => REAL64
  use class_random_uniform, only : uniform
  use class_random_normal , only : normal

  implicit none

  integer , parameter :: npts = 10000
  real(dp), parameter :: xlow = 10.0d0
  real(dp), parameter :: xup  = 20.0d0

  real(dp), allocatable :: xi(:)
  type(uniform) :: urand
  type(normal)  :: nrand
  integer       :: i

  nrand = normal(mu=10.0d0, sigma=1.0d-1)

  allocate(xi(npts))
  call nrand % get(xi)
  print *, "xi ", sum(xi)/dble(npts)

  open(unit=10, file="samples.dat") 
  do i = 1, npts
     write(10,*) xi(i)
  end do
  close(10)

  print *, "xmin, xval ", minval(xi), maxval(xi)
  deallocate(xi)

end subroutine test
