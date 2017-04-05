#include "scalar.fpp"

!=====================================================================!
! This is a simple example of setting up a physical system for
! analysis within the framework.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module test_ode_class

  use constants                 , only : WP
  use dynamic_physics_interface , only : dynamics

  implicit none

  private

  public :: ODE

  !-------------------------------------------------------------------!
  ! Type that implements a test ODE of supplied order
  !-------------------------------------------------------------------!
  
  type, extends(dynamics) :: ODE
     
     ! Define constants and other parameters needed for residual and
     ! jacobian assembly here
     
     type(scalar), allocatable :: A(:)
     
   contains

     ! Implement deferred procedures from superclasses
     procedure :: add_residual
     procedure :: add_jacobian
     procedure :: get_initial_condition

     ! Destructor
     final :: destruct_ODE

  end type ODE

  ! Interface to construct the ODE
  interface ODE
     procedure construct_ODE
  end interface ODE

contains
 
  !===================================================================!
  ! Constructor for ODE system
  !===================================================================!
  
  pure type(ODE) function construct_ODE(A, order, nvars) &
       & result (this)

    type(scalar)  , intent(in) :: A(:)
    type(integer) , intent(in) :: order, nvars

    ! Set the number of state variables of ODE system
    call this % set_num_state_vars(nvars)
    
    ! Set time order of ODE system
    call this % set_time_deriv_order(order)

    ! Set the system parameters
    allocate(this % A, source=A)
   
  end function construct_ODE
  
  !===================================================================!
  ! Destructor for ODE system
  !===================================================================!
  
  pure subroutine destruct_ODE(this)

    type(ODE), intent(inout) :: this

    deallocate(this % A)
    
  end subroutine destruct_ODE
  
  !===================================================================!
  ! Residual assembly at each time step. This is a mandary function
  ! that the user needs to implement to use the integration
  ! scheme. This is where the differential equations are supplied to
  ! the solver.
  ! ===================================================================!
  
  pure subroutine add_residual(this, residual, U)

    class(ODE)   , intent(inout) :: this
    type(scalar) , intent(inout) :: residual(:)
    type(scalar) , intent(in)    :: U(:,:)
    
    assemble: block

      type(integer) :: order

      order = this % get_time_deriv_order()

      residual(:) = residual(:) + U(order+1,:) - (this % A(:)**order) * U(1,:) 

    end block assemble
    
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

    class(ODE) , intent(inout) :: this
    type(scalar) , intent(inout) :: jacobian(:,:)
    type(scalar) , intent(in)    :: coeff(:)
    type(scalar) , intent(in)    :: U(:,:)
    
    assemble: block

      type(integer) :: j, nvars, order
      
      nvars = this % get_num_state_vars()
      order = this % get_time_deriv_order() 
      
      forall(j = 1:nvars)
         jacobian(j,j)=jacobian(j,j) - (this%A(j)**order)*coeff(1) + coeff(order+1)
      end forall

    end block assemble

  end subroutine add_jacobian
  
  !===================================================================!
  ! Sets the initial condition for use in the integator. If first order
  ! system just set initial Q, if a second order system set initial Q
  ! and qdot
  !===================================================================!  

  pure subroutine get_initial_condition(this, U)
    
    class(ODE)   , intent(in)    :: this
    type(scalar) , intent(inout) :: U(:,:)
    type(integer) :: n

    ! Mock the initial conditions 
    forall(n=1:this%time_deriv_order-1) 
       U(n,:) = this % A
    end forall

  end subroutine get_initial_condition

end module test_ode_class
