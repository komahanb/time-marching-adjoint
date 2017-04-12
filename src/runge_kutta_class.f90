#include "scalar.fpp"

!=====================================================================!
! Diagonally Implicit Runge Kutta Integration Module
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================! 

module runge_kutta_integrator_class

  use integrator_interface      , only : integrator
  use dynamic_physics_interface , only : dynamics
  use iso_fortran_env           , only : dp => REAL64

  implicit none

  private
  public :: DIRK
  
  !===================================================================! 
  ! DIRK Integrator type
  !===================================================================! 

  type, extends(integrator) :: DIRK

     ! DIRK variables
     type(integer)             :: max_dirk_order = 4

     ! DIRK Coefficients
     type(scalar), allocatable :: A(:,:)
     type(scalar), allocatable :: B(:)
     type(scalar), allocatable :: C(:)

     ! Tracking variables
     type(integer) :: current_step
     type(integer) :: current_stage

   contains
     
     ! Override integration
     procedure :: integrate => integrate

     ! Helper routines
     procedure :: evaluate_states
     procedure :: get_linear_coeff
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
       & max_dirk_order) result(this)

    class(dynamics)   , intent(in)   , target :: system
    type(scalar)      , intent(in)            :: tinit, tfinal
    type(scalar)      , intent(in)            :: h
    type(integer)     , intent(in)            :: max_dirk_order
    type(logical)     , intent(in)            :: implicit   
    type(integer) :: num_stages

    print *, "======================================"
    print *, ">> Diagonally Implicit Runge Kutta << "
    print *, "======================================"

    !-----------------------------------------------------------------!
    ! Set the order of integration
    !-----------------------------------------------------------------!

    this % max_dirk_order = max_dirk_order
    print '("  >> Max DIRK Order          : ",i4)', this % max_dirk_order

    num_stages = this % max_dirk_order - 1
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
  ! Time integration logic
  !===================================================================!
  
  impure subroutine integrate( this )    

    class(dirk), intent(inout) :: this
    integer :: k, i

    ! Set states to zero
    this % U     = 0.0d0
    this % time  = 0.0d0

    ! Get the initial condition
    call this % system % get_initial_condition(this % U(1,:,:))

    ! March in time
    time: do k = 2, this % num_steps

       this % current_step = k

       stage: do i = 1, this % num_stages

          this % current_stage = i

          !call this % evaluate_time(this % time(k), this % time(k-1), this % C(i)*this % h)
          
          !call this % evaluate_states(this % time(1:k), this % U(1:k,:,:))

       end do stage

       ! Advance the state to the current step
       ! call this % timeMarch(this % u, this % udot, this % uddot)

    end do time

  end subroutine integrate

  !===================================================================!
  ! Routine that checks if the Butcher Tableau entries are valid for
  ! the chosen number of stages/order
  !===================================================================!
  
  subroutine check_coeffs(this)

    class(DIRK) :: this
    type(integer) :: i

    do i = 1, this  % num_stages
       if (abs(this % C(i) - sum(this % A(i,:))) .gt. 0) then
          print *, "WARNING: sum(A(i,j)) != C(i)", i, this % num_stages
       end if
    end do

    if (abs(sum(this % B) - 1.0d0) .gt. 0) then
       print *, "WARNING: sum(B) != 1", this % num_stages
    end if

  end subroutine check_coeffs

 !===================================================================!
  ! Butcher's tableau for DIRK 
  !===================================================================!

  subroutine setup_coeffs(this)

    class(DIRK) :: this
    type(scalar), parameter :: PI = 3.141592653589793_dp
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
  
  impure subroutine destroy(this)

    type(DIRK), intent(inout) :: this

    ! Parent class call
    call this % destruct()

    ! Deallocate DIRK coefficient
    if(allocated(this % A)) deallocate(this % A)
    if(allocated(this % B)) deallocate(this % B)
    if(allocated(this % C)) deallocate(this % C)

  end subroutine destroy

  !================================================================!
  ! Interface to approximate states using the time marching coeffs
  !================================================================!

  impure subroutine evaluate_states(this, t, u)

    use nonlinear_algebra, only : nonlinear_solve

    class(DIRK)   , intent(in)    :: this
    type(scalar) , intent(in)    :: t(:)      ! array of time values
    type(scalar) , intent(inout) :: u(:,:,:)  ! previous values of state variables

    type(scalar)  , allocatable :: lincoeff(:)  ! order of equation + 1
    type(integer) :: k , i, n
    type(scalar)  :: scale
!!$
!!$    ! Pull out the number of time steps of states provided and add one
!!$    ! to point to the current time step
!!$    k = size(u(:,1,1))
!!$
!!$    if (.not. stage_solve) then
!!$
!!$       ! March q to next time step
!!$       forall(m=1:this%nsvars)
!!$          q(k,m) = q(k-1,m) + this % h*sum(this % B(1:this%num_stages) &
!!$               &* this % QDOT(k, 1:this%num_stages, m))
!!$       end forall
!!$
!!$       ! March qdot
!!$       forall(m=1:this%nsvars)
!!$          qdot(k,m) = qdot(k-1,m) + this % h*sum(this % B(1:this%num_stages) &
!!$               &* this % QDDOT(k, 1:this%num_stages, m))
!!$       end forall
!!$
!!$       ! March qddot
!!$       forall(m=1:this%nsvars)
!!$          qddot(k,m) = sum(this % B(1:this%num_stages) &
!!$               &* this % QDDOT(k,1:this%num_stages,m))
!!$       end forall
!!$
!!$    else
!!$
!!$       associate( &
!!$            & p => this % get_order(k), &
!!$            & A => this % A(this % get_order(k),:), &
!!$            & h => t(k) - t(k-1), &
!!$            & torder => this % system % get_time_deriv_order())
!!$
!!$         ! Assume a value for highest order state
!!$       u(k,torder+1,:) = 0
!!$
!!$       ! Find the lower order states based on DIRK formula
!!$       do n = torder, 1, -1
!!$          u(k,n,:) = u(k-1,n,:)
!!$          do i = 1, this % num_stages
!!$             scale = (t(k-i)-t(k-i-1))*A(i+1)
!!$             u(k,n,:) = u(k,n,:) + scale*u(k-i,n+1,:)
!!$          end do
!!$       end do
!!$
!!$       ! guess qddot
!!$       if (i .eq. 1) then ! copy previous global state
!!$          this % QDDOT(k,i,:) = this % UDDOT(k-1,:)
!!$       else ! copy previous local state
!!$          this % QDDOT(k,i,:) = this % QDDOT(k,i-1,:)
!!$       end if
!!$
!!$       ! compute the stage velocity states for the guessed QDDOT
!!$       forall(m = 1 : this % nsvars)
!!$          this % QDOT(k,i,m) = this % udot(k-1,m) &
!!$               & + this % h*sum(this % A(i,:)&
!!$               & * this % QDDOT(k,:, m))
!!$       end forall
!!$
!!$       ! compute the stage states for the guessed QDDOT
!!$       forall(m = 1 : this % nsvars)
!!$          this % Q(k,i,m) = this % u(k-1,m) &
!!$               & + this % h*sum(this % A(i,:)*this % QDOT(k,:, m))
!!$       end forall
!!$
!!$
!!$       ! Perform a nonlinear solution if this is a implicit method
!!$       if ( this % is_implicit() ) then
!!$          allocate(lincoeff(torder + 1))
!!$          call this % get_linear_coeff(lincoeff, p, h)
!!$          call nonlinear_solve(this % system, lincoeff, t(k), u(k,:,:))
!!$          deallocate(lincoeff)
!!$       end if
!!$
!!$     end associate
!!$
!!$  end if

end subroutine evaluate_states

!================================================================!
! Retrieve the coefficients for linearizing the jacobian
!================================================================!

impure subroutine get_linear_coeff(this, lincoeff, stage, h)

  class(DIRK)   , intent(in)    :: this
  type(scalar)  , intent(in)    :: h
  type(scalar)  , intent(inout) :: lincoeff(:)
  type(integer) , intent(in)    :: stage

  type(scalar)  :: a
  type(integer) :: deriv_order,p

  a           = this % A(stage, stage)
  deriv_order = this % system % get_time_deriv_order()

  forall(p = 0:deriv_order)
     lincoeff(p+1) = (a*h)**(deriv_order-p)
  end forall

end subroutine get_linear_coeff

end module runge_kutta_integrator_class
