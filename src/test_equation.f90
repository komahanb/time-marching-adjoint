#include "scalar.fpp"

!=====================================================================!
! This is a simple example of setting up a physical system for time
! dependent analysis.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module test_equation_class

  use constants                 , only : WP
  use dynamic_physics_interface , only : dynamics

  implicit none

  private

  public :: TODE

  !-------------------------------------------------------------------!
  ! Type that implements a test ODE of supplied differential order
  !-------------------------------------------------------------------!
  
  type, extends(dynamics) :: TODE
     
     ! Define constants and other parameters needed for residual and
     ! jacobian assembly here    
     
   contains

     ! Implement deferred procedures from superclasses
     procedure :: add_residual
     procedure :: add_jacobian
     procedure :: get_initial_condition

     ! Destructor
     final :: destruct_TODE

  end type TODE

  ! Interface to construct the ODE
  interface TODE
     procedure construct_TODE
  end interface TODE

contains
 
  !===================================================================!
  ! Constructor for TODE system
  !===================================================================!
  
  pure type(TODE) function construct_TODE() &
       & result (this)

    ! Set a description
    call this % set_description('TEST-ODE')
    
    ! Set the number of state variables of ODE system
    call this % set_num_state_vars(1)
    
    ! Set time order of ODE system
    call this % set_differential_order(1)

  end function construct_TODE
  
  !===================================================================!
  ! Destructor for ODE system
  !===================================================================!
  
  pure subroutine destruct_TODE(this)

    type(TODE), intent(inout) :: this

  end subroutine destruct_TODE
  
  !===================================================================!
  ! Residual assembly at each time step. This is a mandary function
  ! that the user needs to implement to use the integration
  ! scheme. This is where the differential equations are supplied to
  ! the solver.
  ! ===================================================================!
  
  pure subroutine add_residual(this, residual, U, t)

    class(TODE)  , intent(inout) :: this
    type(scalar) , intent(inout) :: residual(:)
    type(scalar) , intent(in)    :: U(:,:)
    type(scalar) , intent(in)    :: t
    
    residual(:) = residual(:) + U(2,:) - cos(t) * U(1,:) 

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

    class(TODE)  , intent(inout) :: this
    type(scalar) , intent(inout) :: jacobian(:,:)
    type(scalar) , intent(in)    :: coeff(:)
    type(scalar) , intent(in)    :: U(:,:)
    type(scalar) , intent(in)    :: t
        
    jacobian = jacobian + coeff(2) - cos(t)*coeff(1)

  end subroutine add_jacobian
  
  !===================================================================!
  ! Sets the initial condition for use in the integator. If first order
  ! system just set initial Q, if a second order system set initial Q
  ! and qdot
  !===================================================================!  

  pure subroutine get_initial_condition(this, U)
    
    class(TODE)  , intent(in)    :: this
    type(scalar) , intent(inout) :: U(:,:)

    U(1,:) = 1.0d0
    
  end subroutine get_initial_condition

end module test_equation_class
