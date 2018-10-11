#include "scalar.fpp"

!=====================================================================!
! One dimensional unsteady tranport physics. 
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module unsteady_transport_class

  use constants                 , only : WP
  use dynamic_physics_interface , only : dynamics

  implicit none

  private

  public :: unsteady_transport

  !-------------------------------------------------------------------!
  ! Type that implements first order transport equations
  !-------------------------------------------------------------------!
  
  type, extends(dynamics) :: unsteady_transport

     type(scalar) :: conv_speed
     type(scalar) :: diff_coeff 
     
   contains
     
     ! Implement deferred procedures from superclasses
     procedure :: add_residual
     procedure :: add_jacobian
     procedure :: get_initial_condition
     
     ! Destructor
     final :: destruct
     
  end type unsteady_transport
  
  ! Interface to construct the physical system
  interface unsteady_transport
     procedure construct_unsteady_transport
  end interface unsteady_transport

contains
 
  !===================================================================!
  ! Constructor for unsteady transport physics
  !===================================================================!
  
  type(unsteady_transport) function construct_unsteady_transport( &
       & diffusion_coeff, convective_velocity) &
       & result (this)

    type(scalar), intent(in) :: diffusion_coeff, convective_velocity

    ! System parameters
    this % conv_speed = convective_velocity
    this % diff_coeff = diffusion_coeff
    
    ! Set time order of physical system
    call this % set_differential_order(1)
    
    ! Set the number of state variables based on spatial
    ! discretization of the governing equations
    ! call this % set_num_state_vars(2)
    stop
    
  end function construct_unsteady_transport
  
  !===================================================================!
  ! Destructor for unsteady transport physics
  !===================================================================!
  
  pure subroutine destruct(this)
    
    type(unsteady_transport), intent(inout) :: this
    
  end subroutine destruct
  
  !===================================================================!
  ! Residual assembly at each time step
  !===================================================================!
  
  pure subroutine add_residual(this, residual, U)
    
    class(unsteady_transport), intent(inout) :: this
    type(scalar)             , intent(inout) :: residual(:)
    type(scalar)             , intent(in)    :: U(:,:)
    
    associate(phi=>U(1,:), phidot=> U(2,:), &
         & vel=>this % conv_speed, gamma => this % diff_coeff)
      
!!$      residual(1) = residual(1) + phidot(1) - q(2)
!!$      
!!$      residual(2) = residual(2) + phidot(2) -  mu * (1.0_wp &
!!$           & - q(1) * q(1)) * q(2) + q(1)

    end associate
    
  end subroutine add_residual

  !===================================================================!
  ! Jacobian assembly at each time step. 
  !===================================================================!
  
  pure subroutine add_jacobian(this, jacobian, coeff, U)

    class(unsteady_transport) , intent(inout) :: this
    type(scalar)              , intent(inout) :: jacobian(:,:)
    type(scalar)              , intent(in)    :: coeff(:)
    type(scalar)              , intent(in)    :: U(:,:)
    
    associate(phi=>U(1,:), phidot=> U(2,:), &
         & vel=>this % conv_speed, gamma => this % diff_coeff, &
         & alpha=>coeff(1), beta=>coeff(2))

!!$
!!$      DRDQ: block
!!$
!!$        ! derivative of first equation
!!$
!!$        jacobian(1,1) = jacobian(1,1) + alpha*0.0_WP
!!$        jacobian(1,2) = jacobian(1,2) - alpha*1.0_WP
!!$
!!$        ! derivative of second equation
!!$
!!$        jacobian(2,1) = jacobian(2,1) + mu*alpha*(1.0d0 + 2.0d0*q(1)*q(2))
!!$        jacobian(2,2) = jacobian(2,2) + mu*alpha*(q(1)*q(1) - 1.0_WP)
!!$
!!$      end block DRDQ
!!$
!!$      DRDPHIDOT: block
!!$
!!$        ! derivative of first equation
!!$
!!$        jacobian(1,1) = jacobian(1,1) + beta*1.0_WP
!!$        jacobian(1,2) = jacobian(1,2) + beta*0.0_WP
!!$
!!$        ! derivative of second equation
!!$
!!$        jacobian(2,1) = jacobian(2,1) + mu*beta*0.0_WP
!!$        jacobian(2,2) = jacobian(2,2) + mu*beta*1.0_WP
!!$
!!$      end block DRDPHIDOT

    end associate

  end subroutine add_jacobian
  
  !===================================================================!
  ! Sets the initial condition for use in the integator. 
  !===================================================================!  

  pure subroutine get_initial_condition(this, U)
    
    class(unsteady_transport), intent(in)    :: this
    type(scalar)             , intent(inout) :: U(:,:)
    type(scalar), parameter :: pi = 4.0_wp*atan(1.0_wp)

    associate(phi=>U(1,:), phidot=> U(2,:), &        
         & vel=>this % conv_speed, gamma => this % diff_coeff)
    end associate
    
  end subroutine get_initial_condition

end module unsteady_transport_class
