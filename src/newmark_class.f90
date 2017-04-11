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
     procedure :: get_linear_coeff

     ! Destructor
     final :: destroy
     
  end type newmark

  interface newmark
     module procedure initialize
  end interface newmark

contains

  !===================================================================!
  ! Initialize the newmark datatype and allocate required variables
  !===================================================================!
  
  type(newmark) function initialize(system, tinit, tfinal, h, implicit, &
       & max_order) result(this)

    class(dynamics)   , intent(in)   , target :: system
    type(scalar)      , intent(in)            :: tinit, tfinal
    type(scalar)      , intent(in)            :: h
    type(integer)     , intent(in)            :: max_order
    type(logical)     , intent(in)            :: implicit   

    print *, "======================================"
    print *, ">>       Newmark Beta Gamma        << "
    print *, "======================================"
    
    call this % construct(system, tinit, tfinal, h, implicit)
    
    if ( this % time_deriv_order .ne. 2 ) then
       print *, " Warning: Newmark-Beta-Gamma method works for second" &
            & // " order systems in current form..."
       stop
    end if

  end function initialize

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
  
  impure subroutine evaluate_states(this, unew, uold)

    use nonlinear_algebra, only : nonlinear_solve

    class(newmark) , intent(in)  :: this
    type(scalar)   , intent(in)  :: uold(:,:,:)  ! previous values of state variables
    type(scalar)   , intent(out) :: unew(:,:)    ! approximated value at current step    
    type(scalar)   , allocatable :: lincoeff(:)  ! order of equation + 1
    type(integer) :: k , i, n
    type(scalar)  :: scale

    ! Pull out the number of time steps of states provided and add one
    ! to point to the current time step
    k = size(uold(:,1,1)) + 1

    !-----------------------------------------------------------------!
    ! Assume a UDDOT for the next time step
    !-----------------------------------------------------------------!

    unew(3,:) = 0

    !-----------------------------------------------------------------!
    ! Approximate UDOT using NBG
    !-----------------------------------------------------------------!

    unew(2,:) = uold(k-1,2,:) 

    scale = this % h * (1.0d0 - this % GAMMA)
    unew(2,:) = unew(2,:) + scale*uold(k-1,3,:) 

    scale = this % h * this % GAMMA
    unew(2,:) = unew(2,:) + scale*unew(3,:) 

    !-----------------------------------------------------------------!
    ! Approximate U using NBG
    !-----------------------------------------------------------------!

    unew(1,:) = uold(k-1,1,:) 

    scale = this % h
    unew(1,:) = unew(1,:) + scale*uold(k-1,2,:) 

    scale = this % h * this % h * (1.0d0 - 2.0d0 * this % BETA )/2.0d0
    unew(1,:) = unew(1,:) + scale*uold(k-1,3,:) 

    scale = this % h * this % h * this % BETA 
    unew(1,:) = unew(1,:) + scale*unew(3,:) 
    
    ! Perform a nonlinear solution if this is a implicit method
    if ( this % is_implicit() ) then
       allocate(lincoeff(this % time_deriv_order + 1))
       call this % get_linear_coeff(lincoeff, this % h)       
       call nonlinear_solve(this % system, lincoeff, &
            & this % h, unew, this % approximate_jacobian)
       deallocate(lincoeff)
    end if

  end subroutine evaluate_states

  !================================================================!
  ! Retrieve the coefficients for linearizing the jacobian
  !================================================================!
  
  impure subroutine get_linear_coeff(this, lincoeff, h)

    class(Newmark) , intent(in)    :: this
    type(scalar)   , intent(in)    :: h ! step size
    type(scalar)   , intent(inout) :: lincoeff(:)   ! order of equation + 1   
   
    lincoeff(1) = this % BETA*h*h
    lincoeff(2) = this % GAMMA*h
    lincoeff(3) = 1.0d0

  end subroutine get_linear_coeff

end module newmark_integrator_class
