#include "scalar.fpp"

!=====================================================================!
! VANDERPOL OSCILLATOR physics. This is a simple example of setting up
! a physical system for analysis within the framework.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module vanderpol_system

  use constants                        , only : WP
  use dynamic_physics_interface , only : dynamics

  implicit none

  private

  public :: vanderpol

  !-------------------------------------------------------------------!
  ! Type that implements vanderpol equations in first order form
  !-------------------------------------------------------------------!
  
  type, extends(dynamics) :: vanderpol_first_order

     ! Define constants and other parameters needed for residual and
     ! jacobian assembly here

     type(scalar)              :: m
     type(scalar), allocatable :: q(:), qdot(:)
     
   contains

     ! Implement deferred procedures from superclasses
     procedure :: add_residual
     procedure :: add_jacobian
     procedure :: set_initial_condition

     ! Destructor
     final :: destruct

  end type vanderpol_first_order

  ! Interface to construct vanderpol system
  interface vanderpol
     procedure construct_first_order
  end interface vanderpol

contains
 
  !===================================================================!
  ! Constructor for vanderpol system
  !===================================================================!
  
  pure type(vanderpol_first_order) function construct_first_order(mass) &
       & result (this)

    type(scalar), intent(in) :: mass

    ! Set the number of state variables of the vanderpol system
    call this % set_num_state_vars(2)
    
    ! Set time order of vanderpol system
    call this % set_time_order(1)

    ! Set the oscillator parameter
    this % m = mass

    ! Allocate state variables
    allocate(this % q   ( this % get_num_state_vars() ))
    allocate(this % qdot( this % get_num_state_vars() ))

    ! Fetch the initial conditions
    call this % set_initial_condition()    
    
  end function construct_first_order
  
  !===================================================================!
  ! Destructor for vanderpol system
  !===================================================================!
  
  pure subroutine destruct(this)

    type(vanderpol_first_order), intent(inout) :: this

    if (allocated(this % q   )) deallocate(this % q   )
    if (allocated(this % qdot)) deallocate(this % qdot)    

  end subroutine destruct
  
  !===================================================================!
  ! Residual assembly at each time step. This is a mandary function
  ! that the user needs to implement to use the integration
  ! scheme. This is where the differential equations are supplied to
  ! the solver.
  !===================================================================!
  
  pure subroutine add_residual(this, residual)

    class(vanderpol_first_order), intent(inout) :: this
    type(scalar)                , intent(inout) :: residual(:)

    residual(1) = residual(1) + this % qdot(1) - this % q(2)
    residual(2) = residual(2) + this % qdot(2) - this % m * (1.0_wp &
         & - this % q(1) * this % q(1)) * this % q(2) + this % q(1)
    
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
  
  pure subroutine add_jacobian(this, jacobian, coeff)

    class(vanderpol_first_order) , intent(inout) :: this
    type(scalar)                 , intent(inout) :: jacobian(:,:)
    type(scalar)                 , intent(in)    :: coeff(:)

    !-----------------------------------------------------------------!
    ! Add dR/dQ
    !-----------------------------------------------------------------!

    DRDQ: block

      type(scalar) :: alpha

      alpha = coeff(1)

      ! derivative of first equation

      jacobian(1,1) = jacobian(1,1) + alpha*0.0_WP
      jacobian(1,2) = jacobian(1,2) - alpha*1.0_WP

      ! derivative of second equation

      jacobian(2,1) = jacobian(2,1) + this % m*alpha*(1.0_WP &
           & + 2.0_WP*this % q(1)*this % q(2))
      jacobian(2,2) = jacobian(2,2) + this % m*alpha*(this % q(1)*this % q(1) - 1.0_WP)

    end block DRDQ

    !-----------------------------------------------------------------!
    ! Add dR/dQDOT
    !-----------------------------------------------------------------!

    DRDQDOT: block

      type(scalar) :: beta

      beta = coeff(2)
      
      ! derivative of first equation

      jacobian(1,1) = jacobian(1,1) + beta*1.0_WP
      jacobian(1,2) = jacobian(1,2) + beta*0.0_WP

      ! derivative of second equation

      jacobian(2,1) = jacobian(2,1) + this % m*beta*0.0_WP
      jacobian(2,2) = jacobian(2,2) + this % m*beta*1.0_WP

    end block DRDQDOT

  end subroutine add_jacobian
  
  !===================================================================!
  ! Sets the initial condition for use in the integator. If first order
  ! system just set initial Q, if a second order system set initial Q
  ! and udot
  !===================================================================!  

  pure subroutine set_initial_condition(this)

    class(vanderpol_first_order), intent(inout) :: this

    this % q(1) = 1.0_WP
    this % q(2) = 2.0_WP
      
  end subroutine set_initial_condition

end module vanderpol_system
