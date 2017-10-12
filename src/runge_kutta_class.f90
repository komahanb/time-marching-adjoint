#include "scalar.fpp"

!=====================================================================!
! Diagonally Implicit Runge Kutta Integration Module
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================! 

module runge_kutta_integrator_class

  use integrator_interface      , only : integrator
  use dynamic_physics_interface , only : dynamics
  use constants                 , only : DP, TINY, PI

  implicit none

  private
  public :: DIRK
  
  !===================================================================! 
  ! DIRK Integrator type
  !===================================================================! 

  type, extends(integrator) :: DIRK

     ! DIRK variables
     type(integer)             :: max_order = 4

     ! DIRK Coefficients
     type(scalar), allocatable :: A(:,:)
     type(scalar), allocatable :: B(:)
     type(scalar), allocatable :: C(:)

   contains     

     ! Helper routines
     procedure :: step    
     procedure :: get_linearization_coeff
     procedure :: get_bandwidth
     procedure :: setup_coeffs
     procedure :: check_coeffs

     ! Destructor
     final :: destroy

  end type DIRK

  interface DIRK
     module procedure create
  end interface DIRK

contains

  !===================================================================!
  ! Initialize the DIRK datatype and allocate required variables
  !===================================================================!
  
  type(dirk) function create(system, tinit, tfinal, h, implicit, &
       & accuracy_order) result(this)

    class(dynamics)   , intent(in)   , target :: system
    type(scalar)      , intent(in)            :: tinit, tfinal
    type(scalar)      , intent(in)            :: h
    type(integer)     , intent(in)            :: accuracy_order
    type(logical)     , intent(in)            :: implicit   
    type(integer) :: num_stages

    print *, "======================================"
    print *, ">> Diagonally Implicit Runge Kutta << "
    print *, "======================================"

    this % max_order = accuracy_order
    print '("  >> Max DIRK Order       : ", i4)', this % max_order

    num_stages = this % max_order - 1
    call this % construct(system, tinit, tfinal, h, implicit, num_stages)

    allocate( this % A (num_stages, num_stages) )
    this % A = 0.0d0 
    allocate( this % B (num_stages) )
    this % B = 0.0d0 
    allocate( this % C (num_stages) )
    this % C = 0.0d0 

    call this % setup_coeffs()
    call this % check_coeffs()

  end function create
  
  !===================================================================!
  ! Returns the order of approximation for the given time step k and
  ! degree d
  !===================================================================!
  
  pure type(integer) function get_bandwidth(this, time_index) result(width)

    class(DIRK)  , intent(in) :: this
    type(integer), intent(in) :: time_index

    integer :: stage, step
    
    ! Identify the step and stage index based on the global time index
    stage = mod(time_index-1,this % num_stages+1)
    step = 1 + (time_index-1)/(this % num_stages+1)

    if ( stage .eq. 0 ) then 
       width = this % num_stages + 1 ! all stages + current 
    else 
       width = stage
    end if

  end function get_bandwidth

  !===================================================================!
  ! Routine that checks if the Butcher Tableau entries are valid for
  ! the chosen number of stages/order
  !===================================================================!
  
  subroutine check_coeffs(this)

    class(DIRK) :: this
    type(integer) :: i

    do i = 1, this  % num_stages
       if (abs(this % C(i) - sum(this % A(i,:))) .gt. TINY) then
          print *, "WARNING: sum(A(i,j)) != C(i)", i, this % num_stages
       end if
    end do

    if (abs(sum(this % B) - 1.0d0) .gt. TINY) then
       print *, "WARNING: sum(B) != 1", this % num_stages
    end if

  end subroutine check_coeffs

  !===================================================================!
  ! Butcher's tableau for DIRK 
  !===================================================================!

  subroutine setup_coeffs(this)

    class(DIRK) :: this
    type(scalar), parameter :: tmp  = 1.0_dp/(2.0_dp*dsqrt(3.0_dp))
    type(scalar), parameter :: half = 1.0_dp/2.0_dp
    type(scalar), parameter :: one  = 1.0_dp
    type(scalar), parameter :: alpha = 2.0_dp*cos(PI/18.0_dp)/dsqrt(3.0_dp)

    ! Put the entries into the tableau (ROGER ALEXANDER 1977)
    if (this % num_stages .eq. 1) then 

       ! Implicit mid-point rule (A-stable)

       this % A(1,1)    = half
       this % B(1)      = one
       this % C(1)      = half

!!$       this % order     = 2
       
!!$       ! Implicit Euler (Backward Euler) but first order accurate
!!$       this % A(1,1) = one
!!$       this % B(1)   = one
!!$       this % C(1)   = one
!!$       this % order  = 1
       
    else if (this % num_stages .eq. 2) then

       ! Crouzeix formula (A-stable)

       this % A(1,1)    = half + tmp
       this % A(2,1)    = -one/dsqrt(3.0_dp)
       this % A(2,2)    = this % A(1,1)

       this % B(1)      = half
       this % B(2)      = half

       this % C(1)      = half + tmp
       this % C(2)      = half - tmp

!!$       this % order     = 3

    else if (this % num_stages .eq. 3) then

       ! Crouzeix formula (A-stable)

       this % A(1,1)    = (one+alpha)*half
       this % A(2,1)    = -half*alpha
       this % A(3,1)    =  one + alpha

       this % A(2,2)    = this % A(1,1)
       this % A(3,2)    = -(one + 2.0_dp*alpha)
       this % A(3,3)    = this % A(1,1)

       this % B(1)      = one/(6.0_dp*alpha*alpha)
       this % B(2)      = one - one/(3.0_dp*alpha*alpha)
       this % B(3)      = this % B(1)

       this % C(1)      = (one + alpha)*half
       this % C(2)      = half
       this % C(3)      = (one - alpha)*half

!!$       this % order     = 4

    else if (this % num_stages .eq. 4) then

       stop "Four stage DIRK formula does not exist"

    else
       
       print *, this % num_stages
       stop "DIRK Butcher tableau is not implemented for the requested&
            & order/stages"

    end if

  end subroutine setup_coeffs

  !=================================================================!
  ! Destructor for the DIRK integrator
  !=================================================================!
  
  pure subroutine destroy(this)

    type(DIRK), intent(inout) :: this

    ! Parent class call
    call this % destruct()

    ! Deallocate DIRK coefficient
    if(allocated(this % A)) deallocate(this % A)
    if(allocated(this % B)) deallocate(this % B)
    if(allocated(this % C)) deallocate(this % C)

  end subroutine destroy
     
  !================================================================!
  ! Take a time step using the supplied time step and order of
  ! accuracy
  ! ================================================================!

  impure subroutine step(this, t, u, h, p, ierr)

    use nonlinear_algebra, only : solve

    ! Argument variables
    class(DIRK)  , intent(inout) :: this
    type(scalar) , intent(inout) :: t(:)
    type(scalar) , intent(inout) :: u(:,:,:)
    type(integer), intent(in)    :: p
    type(scalar) , intent(in)    :: h
    type(integer), intent(out)   :: ierr

    ! Local variables
    type(scalar), allocatable :: lincoeff(:)  ! order of equation + 1
    type(integer) :: torder, n, j, idx
    type(scalar)  :: scale
    type(integer) :: stage_num

    ! Determined the step and stage based on the bandwidth
    stage_num = mod(p,this % num_stages+1)

    ! Determine remaining paramters needed
    idx = size(u(:,1,1))
    torder = this % system % get_differential_order()

    if ( stage_num .eq. 0 ) then

       ! Advance the time to next step
       t(idx) = t(1) + h
       
       ! Highest time derivative
       do j = 1, this % num_stages
          u(idx,torder+1,:) = u(idx,torder+1,:) + this % B(j) * u(idx-j,torder+1,:)
       end do

       ! Find the lower order derivatives
       do n = torder, 1, -1
          u(idx,n,:) = u(1,n,:)
          do j = 1, this % num_stages
             scale = h * this % B(j)
             u(idx,n,:) = u(idx,n,:) + scale*u(idx-j,n+1,:)
          end do
       end do
       
    else
       
       ! Find the new independent coordinate (time)
       t(idx) = t(1) + h * this % C(stage_num)

       ! Assume a value for highest order state
       u(idx,torder+1,:) = 0.0d0

       ! Find the lower order states based on DIRK formula
       do n = torder, 1, -1
          u(idx,n,:) = u(1,n,:)
          do j = 1, stage_num
             scale = h * this % A(stage_num,j)
             u(idx,n,:) = u(idx,n,:) + scale*u(idx-j+1,n+1,:)
          end do
       end do

       ! Perform a nonlinear solution if this is a implicit method
       if ( this % implicit ) then
          allocate(lincoeff(torder + 1))
          call this % get_linearization_coeff(stage_num, h, lincoeff)
          call solve(this % system, lincoeff, t(idx), u(idx,:,:))
          deallocate(lincoeff)
       end if

    end if

  end subroutine step

!================================================================!
! Retrieve the coefficients for linearizing the jacobian
!================================================================!

pure subroutine get_linearization_coeff(this, cindex, h, lincoeff)

  class(DIRK)   , intent(in)    :: this
  type(integer) , intent(in)    :: cindex
  type(scalar)  , intent(in)    :: h
  type(scalar)  , intent(inout) :: lincoeff(:)
  type(scalar)  :: a
  type(integer) :: deriv_order, n

  a = this % A(cindex, cindex)
  deriv_order = this % system % get_differential_order()

  forall(n=0:deriv_order)
     lincoeff(n+1) = (a*h)**(deriv_order-n)
  end forall

end subroutine get_linearization_coeff

end module runge_kutta_integrator_class
