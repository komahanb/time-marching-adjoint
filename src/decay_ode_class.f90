#include "scalar.fpp"

!=====================================================================!
! This is a simple example of setting up a physical system for time
! dependent analysis. We implement first order decay ODE.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module decay_ode_class

  use constants                 , only : WP
  use dynamic_physics_interface , only : dynamics

  implicit none

  private

  public :: decay

  !-------------------------------------------------------------------!
  ! Type that implements a decay ordinary differential equation
  !-------------------------------------------------------------------!
  
  type, extends(dynamics) :: decay
     
     ! Define constants and other parameters needed for residual and
     ! jacobian assembly here
     
     type(scalar) :: gamma
     
   contains

     ! Implement deferred procedures from superclasses
     procedure :: add_residual
     procedure :: add_jacobian
     procedure :: get_initial_condition

     ! Destructor
     final :: destruct

  end type decay
  
  ! Interface to construct the decay
  interface decay
     procedure construct
  end interface decay
  
contains
 
  !===================================================================!
  ! Constructor for decay system
  !===================================================================!
  
  pure type(decay) function construct(gamma) &
       & result (this)

    type(scalar)  , intent(in) :: gamma

    ! Set the number of state variables of decay system
    call this % set_num_state_vars(1)
    
    ! Set time order of decay system
    call this % set_differential_order(1)

    ! Set the system parameters
    this % gamma = gamma
    
  end function construct
  
  !===================================================================!
  ! Destructor for decay system
  !===================================================================!
  
  pure subroutine destruct(this)

    type(decay), intent(inout) :: this
   
  end subroutine destruct
  
  !===================================================================!
  ! Residual assembly at each time step. This is a mandary function
  ! that the user needs to implement to use the integration
  ! scheme. This is where the differential equations are supplied to
  ! the solver.
  ! ===================================================================!
  
  pure subroutine add_residual(this, residual, U, X)

    class(decay) , intent(inout) :: this
    type(scalar) , intent(inout) :: residual(:)
    type(scalar) , intent(in)    :: U(:,:)
    type(scalar) , intent(in)    :: X(:,:)
    
    ! R = qdot + gamma*q
    residual(:) = residual(:) + U(2,:) + this % gamma * U(1,:)
        
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

    class(decay) , intent(inout) :: this
    type(scalar) , intent(inout) :: jacobian(:,:)
    type(scalar) , intent(in)    :: coeff(:)
    type(scalar) , intent(in)    :: U(:,:)
    type(scalar) , intent(in)    :: X(:,:)
    
    jacobian(1,1) = jacobian(1,1) + coeff(2) + coeff(1)*this % gamma

  end subroutine add_jacobian
  
  !===================================================================!
  ! Sets the initial condition for use in the integator. If first order
  ! system just set initial Q, if a second order system set initial Q
  ! and qdot
  !===================================================================!  

  pure subroutine get_initial_condition(this, U, X)
    
    class(decay) , intent(in)    :: this
    type(scalar) , intent(inout) :: U(:,:)
    type(scalar) , intent(in)    :: X(:,:)

    U(1,:) = 1.0d0

  end subroutine get_initial_condition

end module decay_ode_class
