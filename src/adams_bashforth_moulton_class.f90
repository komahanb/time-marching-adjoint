#include "scalar.fpp"

!=====================================================================!
! Adams Bashworth Moulton Integration Module for first and second
! order systems with adjoint derivative capabilities.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================! 

module abm_integrator_class

  use integrator_interface      , only : integrator
  use dynamic_physics_interface , only : dynamics

  implicit none

  private
  public :: ABM
  
  !===================================================================! 
  ! ABM Integrator type
  !===================================================================! 

  type, extends(integrator) :: ABM
     
     private

     ! ABM variables
     type(integer)             :: max_abm_order = 3
     type(scalar), allocatable :: A(:,:)

   contains
           
     procedure :: evaluate_states
     procedure :: get_linear_coeff
     procedure :: get_order

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

    call this % construct(system, tinit, tfinal, h, implicit, 0)

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
  
  impure subroutine destroy(this)

    type(ABM), intent(inout) :: this

    ! Parent class call
    call this % destruct()

    ! Deallocate ABM coefficient
    if(allocated(this % A)) deallocate(this % A)

  end subroutine destroy

  !===================================================================!
  ! Returns the order of approximation for the given time step k and
  ! degree d
  !===================================================================!

  impure type(integer) function get_order(this, step) result(order)

    class(ABM)   , intent(in) :: this
    type(integer), intent(in) :: step

    order = step - 1

    if (order .gt. this % max_abm_order) order = this % max_abm_order

  end function get_order

  !================================================================!
  ! Interface to approximate states using the time marching coeffs
  !================================================================!

  impure subroutine evaluate_states(this, t, u)

    use nonlinear_algebra, only : nonlinear_solve

    class(ABM)   , intent(in)    :: this
    type(scalar) , intent(in)    :: t(:)      ! array of time values
    type(scalar) , intent(inout) :: u(:,:,:)  ! previous values of state variables

    type(scalar)  , allocatable :: lincoeff(:)  ! order of equation + 1
    type(integer) :: k , i, n
    type(scalar)  :: scale
    
    ! Pull out the number of time steps of states provided and add one
    ! to point to the current time step
    k = size(u(:,1,1))

    associate( &
         & p => this % get_order(k), &
         & A => this % A(this % get_order(k),:), &
         & h => t(k) - t(k-1), &
         & torder => this % system % get_time_deriv_order())

      ! Assume a value for highest order state
      u(k,torder+1,:) = 0

      ! Find the lower order states based on ABM formula
      do n = torder, 1, -1
         u(k,n,:) = u(k-1,n,:)! + h*A(1)*u(k,n+1,:)
         do i = 1, p
            scale = h*A(i)
            u(k,n,:) = u(k,n,:) + scale*u(k+1-i,n+1,:)
         end do
      end do

      ! Perform a nonlinear solution if this is a implicit method
      if ( this % is_implicit() ) then
         allocate(lincoeff(torder + 1))
         call this % get_linear_coeff(lincoeff, p, h)           
         call nonlinear_solve(this % system, lincoeff, t(k), u(k,:,:))
         deallocate(lincoeff)
      end if

    end associate

end subroutine evaluate_states

  !================================================================!
  ! Retrieve the coefficients for linearizing the jacobian
  !================================================================!
  
  impure subroutine get_linear_coeff(this, lincoeff, int_order, h)

    class(ABM)    , intent(in)  :: this
    type(integer) , intent(in)  :: int_order     ! order of approximation of the integration
    type(scalar)  , intent(in)  :: h ! step size
    type(scalar)  , intent(inout) :: lincoeff(:)   ! order of equation + 1   
    type(integer) :: p
    
    associate(&
         & deriv_order => this % system % get_time_deriv_order(), &
         & a => this % A(int_order,1))
      
      forall(p = 0:deriv_order)
         lincoeff(p+1) = (a*h)**(deriv_order-p)
      end forall
      
    end associate

  end subroutine get_linear_coeff

end module abm_integrator_class
