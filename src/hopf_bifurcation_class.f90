#include "scalar.fpp"

!=====================================================================!
! Hopf Bifurcation physics. This is a simple example of setting up
! a physical system for analysis within the framework.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module hopf_bifurcation_system

  use constants                 , only : WP
  use dynamic_physics_interface , only : dynamics

  implicit none

  private

  public :: hopf

  !-------------------------------------------------------------------!
  ! Type that implements hopf bifurcation equations in first order form
  !-------------------------------------------------------------------!
  
  type, extends(dynamics) :: hopf

     ! Define constants and other parameters needed for residual and
     ! jacobian assembly here

     type(scalar) :: a, b
     
   contains

     ! Implement deferred procedures from superclasses
     procedure :: add_residual
     procedure :: add_jacobian
     procedure :: get_initial_condition

     ! Destructor
     final :: destruct

  end type hopf

  ! Interface to construct birfurcation system
  interface hopf
     procedure construct_hopf
  end interface hopf

contains
 
  !===================================================================!
  ! Constructor for Hopf system
  !===================================================================!
  
  pure type(hopf) function construct_hopf(a, b) &
       & result (this)

    type(scalar), intent(in) :: a, b

    ! Set the number of state variables of the bifurcation system
    call this % set_num_state_vars(2)
    
    ! Set time order of bifurcation system
    call this % set_differential_order(1)

    ! Set the oscillator parameter
    this % a = a
    this % b = b
    
  end function construct_hopf
  
  !===================================================================!
  ! Destructor for Hopf system
  !===================================================================!
  
  pure subroutine destruct(this)

    type(hopf), intent(inout) :: this
    
  end subroutine destruct
  
  !===================================================================!
  ! Residual assembly at each time step. This is a mandary function
  ! that the user needs to implement to use the integration
  ! scheme. This is where the differential equations are supplied to
  ! the solver.
  !===================================================================!
  
  pure subroutine add_residual(this, residual, U, t)

    class(hopf)  , intent(inout) :: this
    type(scalar) , intent(inout) :: residual(:)
    type(scalar) , intent(in)    :: U(:,:)
    type(scalar) , intent(in)    :: t

    associate(y => U(1,:), ydot => U(2,:), a => this%a, b => this%b)

      residual(1) = residual(1) + ydot(1) - a + y(1) + 4.0d0*y(1)*y(2)/(1.0d0 + y(1)*y(1))
      residual(2) = residual(2) + ydot(2) - b*y(1)*(1.0d0-y(2)/(1.0d0 + y(1)*y(1)))

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
  
  pure subroutine add_jacobian(this, jacobian, coeff, U, t)

    class(hopf)  , intent(inout) :: this
    type(scalar) , intent(inout) :: jacobian(:,:)
    type(scalar) , intent(in)    :: coeff(:)
    type(scalar) , intent(in)    :: U(:,:)
    type(scalar) , intent(in)    :: t
    
    associate(y => U(1,:), ydot => U(2,:), &
         & a => this % a, b => this % b, &
         & alpha => coeff(1), beta => coeff(2) &
         & )

      DRDQ: block

        ! derivative of first equation
               
        jacobian(1,1) = jacobian(1,1) + alpha*((y(1)**4 + (2.0d0 - 4.0d0*y(2))*y(1)**2 + 1.0d0)/(1.0d0 + y(1)**2)**2)
        jacobian(1,2) = jacobian(1,2) + alpha*(4.0d0*y(1)/(1.0d0 + y(1)**2))

        ! derivative of second equation
        
        jacobian(2,1) = jacobian(2,1) - alpha*(b*(y(1)**4 + (y(2) + 2.0d0)*y(1)**2 - y(2) + 1.0d0)/(1.0d0 + y(1)**2)**2)
        jacobian(2,2) = jacobian(2,2) + alpha*(b*y(1)/(1.0d0+y(1)**2))        

      end block DRDQ

      DRDQDOT: block

        ! derivative of first equation

        jacobian(1,1) = jacobian(1,1) + beta*1.0_WP
        jacobian(1,2) = jacobian(1,2) + beta*0.0_WP

        ! derivative of second equation

        jacobian(2,1) = jacobian(2,1) + beta*0.0_WP
        jacobian(2,2) = jacobian(2,2) + beta*1.0_WP

      end block DRDQDOT

    end associate

  end subroutine add_jacobian
  
  !===================================================================!
  ! Sets the initial condition for use in the integator. If first order
  ! system just set initial Q, if a second order system set initial Q
  ! and qdot
  !===================================================================!  

  pure subroutine get_initial_condition(this, U)
    
    class(hopf)  , intent(in) :: this
    type(scalar) , intent(inout) :: U(:,:)

    U(1,1:2) = [ 0.0_WP, 2.0_WP ]
    
  end subroutine get_initial_condition

end module hopf_bifurcation_system
