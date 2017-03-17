#include "scalar.fpp"

!=====================================================================!
! Adams Bashworth Moulton Integration Module for first and second
! order systems with adjoint derivative capabilities.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================! 

module abm_integrator_class

  use integrator_interface, only : integrator
  use dynamic_physics_interface, only : dynamics

  implicit none

  private
  public :: ABM
  
  type(scalar), parameter :: ONE  = 1.0d0
  type(scalar), parameter :: ZERO = 0.0d0

  !===================================================================! 
  ! ABM Integrator type
  !===================================================================! 

  type, extends(integrator) :: ABM
     
     private

     ! ABM variables
     type(integer) :: max_abm_order = 3
     type(scalar), allocatable, dimension(:,:) :: A(:,:)

   contains
           
     procedure :: evaluate_states
     procedure :: get_state_coeff

     ! Routines for integration
     procedure, private :: approximateStates
     procedure, private :: getLinearCoeff
     procedure, private :: getOrder

     ! Destructor
     final :: destroy

  end type ABM

  interface ABM
     module procedure initialize
  end interface ABM

contains

  !===================================================================!
  ! Initialize the ABM datatype and allocate required variables
  !===================================================================!
  
  type(abm) function initialize(system, tinit, tfinal, h, implicit, &
       & max_abm_order) result(this)

    class(dynamics)   , intent(in)   , target :: system
    type(scalar)      , intent(in)            :: tinit, tfinal
    type(scalar)      , intent(in)            :: h
    type(integer)     , intent(in)            :: max_abm_order
    type(logical)     , intent(in)            :: implicit   

    print *, "======================================"
    print *, ">>   Adams Bashforth Moulton       << "
    print *, "======================================"

    call this % construct(system, tinit, tfinal, h, implicit)

    !-----------------------------------------------------------------!
    ! Set the order of integration
    !-----------------------------------------------------------------!

    this % max_abm_order = max_abm_order
    print '("  >> Max ABM Order          : ",i4)', this % max_abm_order

    allocate( this % A (this % max_abm_order, this % max_abm_order) )
    this % A = 0.0d0 

    ! Set the coefficients
    if ( this % max_abm_order .eq. 1 ) then       
       this % A(1,1:1) = (/ 1.0d0 /)
    else if ( this % max_abm_order .eq. 2 ) then
       this % A(1,1:1) = (/ 1.0d0 /)
       this % A(2,1:2) = (/ 1.0d0/2.0d0, 1.0d0/2.0d0 /)
    else if ( this % max_abm_order .eq. 3 ) then
       this % A(1,1:1) = (/ 1.0d0 /)
       this % A(2,1:2) = (/ 1.0d0/2.0d0, 1.0d0/2.0d0 /)
       this % A(3,1:3) = (/ 5.0d0/12.0d0, 8.0d0/12.0d0, -1.0d0/12.0d0 /)
    else 
       print *,  "Wrong max_abm_order:", this % max_abm_order
       stop
    end if

    ! Sanity check on ABM coeffs
    sanity_check: block
      type(integer) :: j
      do j = 1, this % max_abm_order
         if ( real(sum(this % A(j,1:j)) - 1.0d0) .gt. 0.00001 ) then
            stop "Error in ABM Coeff"
         end if
      end do
    end block sanity_check

  end function initialize

  !=================================================================!
  ! Destructor for the ABM integrator
  !=================================================================!
  
  pure subroutine destroy(this)

    type(ABM), intent(inout) :: this

    ! Parent class call
    call this % destruct()

    ! Deallocate ABM coefficient
    if(allocated(this % A)) deallocate(this % A)

  end subroutine destroy

  !===================================================================!
  ! Returns the linearization scalar coefficients: alpha, beta, gamma
  !===================================================================!

  subroutine getLinearCoeff( this, alpha, beta, gamma )

    class(ABM)   :: this
    type(scalar), intent(out) :: alpha, beta, gamma
    integer :: k, m
!!$
!!$    k = this % current_step
!!$    m = this % getOrder(k)
!!$
!!$    if ( this % second_order ) then
!!$       gamma = 1.0d0
!!$       beta  = this % A(m,1) * this % h
!!$       alpha = beta * this % A(m,1) * this % h
!!$    else 
!!$       gamma = 0.0d0
!!$       beta  = 1.0d0
!!$       alpha = this  % A(m,1) * this % h
!!$    end if

  end subroutine getLinearCoeff

  !===================================================================!
  ! Returns the order of approximation for the given time step k and
  ! degree d
  !===================================================================!

  pure integer function getOrder(this, k)

    class(ABM), intent(in) :: this
    integer, intent(in)    :: k
!!$
!!$    ! Find the order of approximation
!!$    getOrder = k - 1
!!$
!!$    ! Do not let exceed the max order sought
!!$    if ( getOrder .gt. this % max_abm_order ) getOrder = this % max_abm_order

  end function getOrder

  !===================================================================!
  ! Approximate the state variables at each step using ABM formulae
  !===================================================================!

  subroutine approximateStates( this )

    class(ABM)   :: this
    integer      :: k, m, i
    type(scalar) :: scale
!!$
!!$    k = this % current_step
!!$
!!$    m = this % getOrder(k)
!!$
!!$    ! Approximate UDDOT
!!$    this % uddot(k,:) = this % uddot(k-1,:)
!!$
!!$    ! Approximate UDOT
!!$    this % udot(k,:) = this % udot(k-1,:)
!!$
!!$    do i = 0, m-1
!!$       scale = this % h * this % A(m,i+1)
!!$       this % udot(k,:) = this % udot(k,:) + scale * this % uddot(k-i,:)
!!$    end do
!!$
!!$    ! Approximate U
!!$    this % u(k,:) = this % u(k-1,:)
!!$
!!$    do i = 0, m-1
!!$       scale = this % h * this % A(m,i+1)
!!$       this % u(k,:) = this % u(k,:) + scale * this % udot(k-i,:)
!!$    end do

  end subroutine approximateStates

    !================================================================!
  ! Interface to approximate states using the time marching coeffs
  !================================================================!
  
  pure subroutine evaluate_states(this, unew, uold)

    class(ABM), intent(in)    :: this
    type(scalar)     , intent(in)    :: uold(:,:,:)  ! previous values of state variables
    type(scalar)     , intent(inout)   :: unew(:,:)    ! approximated value at current step
    
    type(integer) :: k, i
    type(integer) :: m
    type(scalar) :: scale
     
    ! Pull out the number of time steps of states provided and add one
    ! to point to the current time step
    k = size(uold(:,1,1)) + 1
    
    associate(q=>uold(1,:,:),qdot=>uold(2,:,:),u=>unew(1,:),udot=>unew(2,:))
      
      ! Get the order based on the index 'k'
      ! m = this % get_order(k)
      
      ! Approximate UDOT
      udot = qdot(k-1,:)
      
      ! Approximate U
      u = q(k-1,:)
      do i = 0, m-1
         scale = this % h * this % A(m,i+1)
         u = u + scale * qdot(k-i,:)
      end do

    end associate

  end subroutine evaluate_states

  !================================================================!
  ! Retrieve the state approximation coefficients
  !================================================================!
  
  pure subroutine get_state_coeff(this, scoeff, int_order, time_deriv_order)

    class(ABM) , intent(in)    :: this
    type(integer)     , intent(in)    :: int_order(:)     ! order of approximation of the integration
    type(integer)     , intent(in)    :: time_deriv_order !  order of the differential equation in time
    type(scalar)      , intent(out)   :: scoeff(:)        ! order of equation + 1

    ! scoeff(1:time_deriv_order+1) = this % A(int_order,1:time_deriv_order+1)

  end subroutine get_state_coeff

end module abm_integrator_class
