#include "scalar.fpp"

!=====================================================================!
! This is a simple example of setting up a physical system for time
! dependent analysis.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module first_order_ode_class

  use constants                 , only : WP
  use dynamic_physics_interface , only : dynamics

  implicit none

  private

  public :: FODE

  !-------------------------------------------------------------------!
  ! Type that implements a test FODE of supplied differential order
  !-------------------------------------------------------------------!
  
  type, extends(dynamics) :: FODE
     
     ! Define constants and other parameters needed for residual and
     ! jacobian assembly here
     
     type(scalar) :: a
     type(scalar) :: b     
     
   contains

     ! Implement deferred procedures from superclasses
     procedure :: add_residual
     procedure :: add_jacobian
     procedure :: get_initial_condition

     ! Destructor
     final :: destruct_FODE

  end type FODE

  ! Interface to construct the FODE
  interface FODE
     procedure construct_FODE
  end interface FODE

contains
 
  !===================================================================!
  ! Constructor for FODE system
  !===================================================================!
  
  pure type(FODE) function construct_FODE(a, b) &
       & result (this)

    type(scalar)  , intent(in) :: a, b

    ! Set the number of state variables of FODE system
    call this % set_num_state_vars(1)

    ! Set time order of FODE system
    call this % set_differential_order(1)

    ! Set the system parameters
    this % a = a
    this % b = b

  end function construct_FODE
  
  !===================================================================!
  ! Destructor for FODE system
  !===================================================================!
  
  pure subroutine destruct_FODE(this)

    type(FODE), intent(inout) :: this

  end subroutine destruct_FODE
  
  !===================================================================!
  ! Residual assembly at each time step. This is a mandary function
  ! that the user needs to implement to use the integration
  ! scheme. This is where the differential equations are supplied to
  ! the solver.
  ! ===================================================================!
  
  pure subroutine add_residual(this, residual, U, X)

    class(FODE)  , intent(inout) :: this
    type(scalar) , intent(inout) :: residual(:)
    type(scalar) , intent(in)    :: U(:,:)
    type(scalar) , intent(in)    :: X(:,:)
    
    residual = residual + U(2,:) - this % b * U(1,:) - this % a

  end subroutine add_residual

  !===================================================================!
  ! Jacobian assembly at each time step. 
  !
  ! Jacobian is the matrix of partial derivatives. Each row in the
  ! Jacobian matrix arises from differentiating a single
  ! equation. Each column in the Jacobian comes from a variable in the
  ! problem. Note the NEQN should be equal to NVARS for the system to
  ! be solved.
  !
  ! Note: alpha, beta and gamma are scalars that need to be multiplied
  ! with the partial derivatives DRDQ, DRDudot and DRDuddot
  ! respectively.
  !===================================================================!
  
  pure subroutine add_jacobian(this, jacobian, coeff, U, X)

    class(FODE)   , intent(inout) :: this
    type(scalar) , intent(inout) :: jacobian(:,:)
    type(scalar) , intent(in)    :: coeff(:)
    type(scalar) , intent(in)    :: U(:,:)
    type(scalar) , intent(in)    :: X(:,:)

    jacobian(1,1) = jacobian(1,1) + coeff(2) - coeff(1)*this % b

  end subroutine add_jacobian
  
  !===================================================================!
  ! Sets the initial condition for use in the integator. If first order
  ! system just set initial Q, if a second order system set initial Q
  ! and qdot
  !===================================================================!  

  pure subroutine get_initial_condition(this, U, X)
    
    class(FODE)   , intent(in)    :: this
    type(scalar) , intent(inout) :: U(:,:)
    type(scalar) , intent(in)    :: X(:,:)
        
    U(1,:) = 48.0d0

  end subroutine get_initial_condition

end module first_order_ode_class
