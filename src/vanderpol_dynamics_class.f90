#include "scalar.fpp"

!=====================================================================!
! VANDERPOL OSCILLATOR physics. This is a simple example of setting up
! a physical system for analysis within the framework.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module vanderpol_system

  use constants                 , only : WP
  use dynamic_physics_interface , only : dynamics

  implicit none

  private

  public :: vanderpol_first_order

  !-------------------------------------------------------------------!
  ! Type that implements vanderpol equations in first order form
  !-------------------------------------------------------------------!
  
  type, extends(dynamics) :: vanderpol_first_order

     ! Define constants and other parameters needed for residual and
     ! jacobian assembly here

     type(scalar) :: mu
     
   contains

     ! Implement deferred procedures from superclasses
     procedure :: add_residual
     procedure :: add_jacobian
     procedure :: get_initial_condition

     ! Destructor
     final :: destruct

  end type vanderpol_first_order

  ! Interface to construct vanderpol system
  interface vanderpol_first_order
     procedure construct_first_order
  end interface vanderpol_first_order

contains
 
  !===================================================================!
  ! Constructor for vanderpol system
  !===================================================================!
  
  pure type(vanderpol_first_order) function construct_first_order(mu) &
       & result (this)

    type(scalar), intent(in) :: mu

    ! Set the number of state variables of the vanderpol system
    call this % set_num_state_vars(2)
    
    ! Set time order of vanderpol system
    call this % set_time_deriv_order(1)

    ! Set the oscillator parameter
    this % mu = mu
    
  end function construct_first_order
  
  !===================================================================!
  ! Destructor for vanderpol system
  !===================================================================!
  
  pure subroutine destruct(this)

    type(vanderpol_first_order), intent(inout) :: this
    
  end subroutine destruct
  
  !===================================================================!
  ! Residual assembly at each time step. This is a mandary function
  ! that the user needs to implement to use the integration
  ! scheme. This is where the differential equations are supplied to
  ! the solver.
  !===================================================================!
  
  pure subroutine add_residual(this, residual, U)

    class(vanderpol_first_order), intent(inout) :: this
    type(scalar)                , intent(inout) :: residual(:)
    type(scalar)                , intent(in)    :: U(:,:)

    associate( q => U(1,:), qdot => U(2,:), mu => this%mu)

      residual(1) = residual(1) + qdot(1) - q(2)

      residual(2) = residual(2) + qdot(2) -  mu * (1.0_wp &
           & - q(1) * q(1)) * q(2) + q(1)

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
  
  pure subroutine add_jacobian(this, jacobian, coeff, U)

    class(vanderpol_first_order) , intent(inout) :: this
    type(scalar)                 , intent(inout) :: jacobian(:,:)
    type(scalar)                 , intent(in)    :: coeff(:)
    type(scalar)                 , intent(in)    :: U(:,:)
    
    associate(q=>U(1,:), qdot=> U(2,:), mu=>this%mu, alpha=>coeff(1), beta=>coeff(2))

      DRDQ: block

        ! derivative of first equation

        jacobian(1,1) = jacobian(1,1) + alpha*0.0_WP
        jacobian(1,2) = jacobian(1,2) - alpha*1.0_WP

        ! derivative of second equation

        jacobian(2,1) = jacobian(2,1) + mu*alpha*(1.0d0 + 2.0d0*q(1)*q(2))
        jacobian(2,2) = jacobian(2,2) + mu*alpha*(q(1)*q(1) - 1.0_WP)

      end block DRDQ

      DRDQDOT: block

        ! derivative of first equation

        jacobian(1,1) = jacobian(1,1) + beta*1.0_WP
        jacobian(1,2) = jacobian(1,2) + beta*0.0_WP

        ! derivative of second equation

        jacobian(2,1) = jacobian(2,1) + mu*beta*0.0_WP
        jacobian(2,2) = jacobian(2,2) + mu*beta*1.0_WP

      end block DRDQDOT

    end associate

  end subroutine add_jacobian
  
  !===================================================================!
  ! Sets the initial condition for use in the integator. If first order
  ! system just set initial Q, if a second order system set initial Q
  ! and qdot
  !===================================================================!  

  pure subroutine get_initial_condition(this, U)
    
    class(vanderpol_first_order), intent(in)    :: this
    type(scalar)                , intent(inout) :: U(:,:)

    U(1,1:2) = [ 2.0_WP, 0.0_WP ]
    
  end subroutine get_initial_condition

end module vanderpol_system
