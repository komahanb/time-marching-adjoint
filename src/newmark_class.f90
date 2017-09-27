#include "scalar.fpp"

!=====================================================================!
! Newmark Beta Gamma Integration Module for first and second
! order systems.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================! 

module newmark_integrator_class

  use integrator_interface      , only : integrator
  use dynamic_physics_interface , only : dynamics

  implicit none

  private
  public :: newmark
  
  !===================================================================! 
  ! Newmark Integrator type
  !===================================================================! 

  type, extends(integrator) :: newmark
     
     ! Average Constant Accelearation (second order unconditionally stable)
     type(scalar) :: BETA  = 0.25d0
     type(scalar) :: GAMMA = 0.50d0

   contains
           
     procedure :: evaluate_states
     procedure :: get_linearization_coeff

     ! Destructor
     final :: destroy
     
  end type newmark

  interface newmark
     module procedure create
  end interface newmark

contains

  !===================================================================!
  ! Initialize the newmark datatype and allocate required variables
  !===================================================================!
  
  type(newmark) function create(system, tinit, tfinal, h, implicit, &
       & accuracy_order) result(this)

    class(dynamics)   , intent(in)   , target :: system
    type(scalar)      , intent(in)            :: tinit, tfinal
    type(scalar)      , intent(in)            :: h
    type(integer)     , intent(in)            :: accuracy_order
    type(logical)     , intent(in)            :: implicit   

    print *, "======================================"
    print *, ">>       Newmark Beta Gamma        << "
    print *, "======================================"
    
    call this % construct(system, tinit, tfinal, h, implicit, 0)
    
    if ( this % system % get_differential_order() .ne. 2 ) then
       print *, " Warning: Newmark-Beta-Gamma method works for second" &
            & // " order systems in current form..."
       stop
    end if

  end function create

  !=================================================================!
  ! Destructor for the Newmark integrator
  !=================================================================!
  
  impure subroutine destroy(this)

    type(newmark), intent(inout) :: this

    ! Parent class call
    call this % destruct()

  end subroutine destroy

  !================================================================!
  ! Interface to approximate states using the time marching coeffs
  !================================================================!
  
  impure subroutine evaluate_states(this, t, u)

    use nonlinear_algebra, only : nonlinear_solve

    class(newmark) , intent(in)    :: this
    type(scalar)   , intent(in)    :: t(:)      ! array of time values
    type(scalar)   , intent(inout) :: u(:,:,:)  ! previous values of state variables

    type(scalar)   , allocatable :: lincoeff(:)  ! order of equation + 1
    type(integer) :: k , i, n
    type(scalar)  :: scale, h

    ! Determine the step number
    k = size(u(:,1,1))

    ! Determine the current step size
    h = t(k) - t(k-1)

    !-----------------------------------------------------------------!
    ! Assume a UDDOT for the next time step
    !-----------------------------------------------------------------!

    u(k,3,:) = 0.0d0

    !-----------------------------------------------------------------!
    ! Approximate UDOT using NBG
    !-----------------------------------------------------------------!

    u(k,2,:) = u(k-1,2,:) 
    
    scale = h*(1.0d0 - this % GAMMA)
    u(k,2,:) = u(k,2,:) + scale*u(k-1,3,:) 

    scale = h*this % GAMMA
    u(k,2,:) = u(k,2,:) + scale*u(k,3,:) 

    !-----------------------------------------------------------------!
    ! Approximate U using NBG
    !-----------------------------------------------------------------!

    u(k,1,:) = u(k-1,1,:) 

    scale = h
    u(k,1,:) = u(k,1,:) + scale*u(k-1,2,:) 

    scale = h*h*(1.0d0 - 2.0d0 * this % BETA )/2.0d0
    u(k,1,:) = u(k,1,:) + scale*u(k-1,3,:) 

    scale = h*h*this % BETA 
    u(k,1,:) = u(k,1,:) + scale*u(k,3,:)
    
    ! Perform a nonlinear solution if this is a implicit method
    if ( this % is_implicit() ) then
       allocate(lincoeff(this % system % get_differential_order() + 1))
       call this % get_linearization_coeff(lincoeff, h)       
       call nonlinear_solve(this % system, lincoeff, this % time(k), u(k,:,:))
       deallocate(lincoeff)
    end if

  end subroutine evaluate_states

  !================================================================!
  ! Retrieve the coefficients for linearizing the jacobian
  !================================================================!
  
  impure subroutine get_linearization_coeff(this, lincoeff, h)

    class(Newmark) , intent(in)    :: this
    type(scalar)   , intent(in)    :: h ! step size
    type(scalar)   , intent(inout) :: lincoeff(:)   ! order of equation + 1   
   
    lincoeff(1) = this % BETA*h*h
    lincoeff(2) = this % GAMMA*h
    lincoeff(3) = 1.0d0

  end subroutine get_linearization_coeff

end module newmark_integrator_class
