#include "scalar.fpp"

!=====================================================================!
! A particle falling freely under gravity physics. This is a simple
! example of setting up a physical system for analysis within the
! framework.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module freefall_dynamics_class

  use constants                 , only : WP
  use dynamic_physics_interface , only : dynamics

  implicit none

  private

  public :: freefall

  !-------------------------------------------------------------------!
  ! Type that implements freefall mass damper ODE
  !-------------------------------------------------------------------!
  
  type, extends(dynamics) :: freefall

     ! Define constants and other parameters needed for residual and
     ! jacobian assembly here

     type(scalar) :: mass
     type(scalar) :: gravity
     
   contains

     ! Implement deferred procedures from superclasses
     procedure :: add_residual
     procedure :: add_jacobian
     procedure :: get_initial_condition

     ! Destructor
     final :: destruct
     
  end type freefall

  ! Interface to construct freefall particle system
  interface freefall
     procedure construct
  end interface freefall

contains
 
  !===================================================================!
  ! Constructor for freefall system
  !===================================================================!
  
  pure type(freefall) function construct(M, G) &
       & result (this)

    type(scalar), intent(in) :: M, G

    ! Set the number of state variables of FREEFALL system
    call this % set_num_state_vars(1)
    
    ! Set time order of FREEFALL system
    call this % set_differential_order(2)

    ! Set the system parameters
    this % mass = M
    this % gravity = G
    
  end function construct
  
  !===================================================================!
  ! Destructor for freefall system
  !===================================================================!
  
  pure subroutine destruct(this)

    type(freefall), intent(inout) :: this
    
  end subroutine destruct
  
  !===================================================================!
  ! Residual assembly at each time step. This is a mandary function
  ! that the user needs to implement to use the integration
  ! scheme. This is where the differential equations are supplied to
  ! the solver.
  !===================================================================!
  
  pure subroutine add_residual(this, residual, U, X)

    class(freefall)   , intent(inout) :: this
    type(scalar) , intent(inout) :: residual(:)
    type(scalar) , intent(in)    :: U(:,:)
    type(scalar) , intent(in)    :: X(:,:)

    associate( q => U(1,:), qdot => U(2,:), qddot => U(3,:), &
         & m => this % mass, g => this % gravity)
      
      residual = residual + m*qddot - m*g

    end associate
    
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

    class(freefall)   , intent(inout) :: this
    type(scalar) , intent(inout) :: jacobian(:,:)
    type(scalar) , intent(in)    :: coeff(:)
    type(scalar) , intent(in)    :: U(:,:)
    type(scalar) , intent(in)    :: X(:,:)
    
    associate(&
         & q=>U(1,:), qdot=> U(2,:), qddot=> U(3,:), &
         & M=>this%mass, &
         & alpha=>coeff(1), beta=>coeff(2), gamma=>coeff(3)&
         & )
      
      jacobian = jacobian + gamma*M
      
    end associate

  end subroutine add_jacobian
  
  !===================================================================!
  ! Sets the initial condition for use in the integator. If first order
  ! system just set initial Q, if a second order system set initial Q
  ! and qdot
  !===================================================================!  

  pure subroutine get_initial_condition(this, U, X)
    
    class(freefall)   , intent(in)    :: this
    type(scalar) , intent(inout) :: U(:,:)
    type(scalar) , intent(in)    :: X(:,:)
    
    U(1,:) = 0.0_WP ! position
    U(2,:) = 1.0_WP ! velocity
    
  end subroutine get_initial_condition

end module freefall_dynamics_class
