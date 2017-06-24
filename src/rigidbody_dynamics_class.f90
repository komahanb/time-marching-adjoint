#include "scalar.fpp"

!=====================================================================!
! A particle falling freely under gravity physics. This is a simple
! example of setting up a physical system for analysis within the
! framework.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module rigidbody_dynamics_class

  use constants                 , only : WP
  use dynamic_physics_interface , only : dynamics

  implicit none

  private

  public :: rigidbody

  !-------------------------------------------------------------------!
  ! Type that implements rigidbody mass damper ODE
  !-------------------------------------------------------------------!
  
  type, extends(dynamics) :: rigidbody

     ! Define constants and other parameters needed for residual and
     ! jacobian assembly here
     
     ! Inertial properties of the rigidbody (assuming 3d space)
     type(scalar) :: mass
     type(scalar) :: first_moment_mass(3)
     type(scalar) :: second_moment_mass(3,3)

     ! Initial conditions of the rigidbody (assuming 3d space)
     type(scalar) :: position(3)
     type(scalar) :: orientation(3)
     type(scalar) :: linear_velocity(3)
     type(scalar) :: angular_velocity(3)

   contains

     ! Implement deferred procedures from superclasses
     procedure :: add_residual
     procedure :: add_jacobian
     procedure :: get_initial_condition

     ! Destructor
     final :: destruct
     
  end type rigidbody

  ! Interface to construct rigidbody particle system
  interface rigidbody
     procedure construct
  end interface rigidbody

contains
 
  !===================================================================!
  ! Constructor for rigidbody system
  !===================================================================!
  
  pure type(rigidbody) function construct(mass, fmass, smass, &
       & r, theta, v, omega) result (this)
    
    type(scalar), intent(in) :: mass, fmass(:), smass(:,:)
    type(scalar), intent(in) :: r(:), theta(:), v(:), omega(:)
    
    ! Set the number of state variables of RIGIDBODY system
    call this % set_num_state_vars(6)
    
    ! Set time order of RIGIDBODY system
    call this % set_time_deriv_order(2)

    ! Set the object attributes
    this % mass = mass
    this % first_moment_mass = fmass
    this % second_moment_mass = smass

    ! Set the states
    this % position = r
    this % orientation = theta
    this % linear_velocity = v
    this % angular_velocity = omega

  end function construct
  
  !===================================================================!
  ! Destructor for rigidbody system
  !===================================================================!
  
  pure subroutine destruct(this)

    type(rigidbody), intent(inout) :: this
    
  end subroutine destruct
  
  !===================================================================!
  ! Residual assembly at each time step. This is a mandary function
  ! that the user needs to implement to use the integration
  ! scheme. This is where the differential equations are supplied to
  ! the solver.
  !===================================================================!
  
  pure subroutine add_residual(this, residual, U)

    class(rigidbody), intent(inout) :: this
    type(scalar), intent(inout) :: residual(:)
    type(scalar), intent(in)    :: U(:,:)

    associate( q => U(1,:), qdot => U(2,:), qddot => U(3,:), &
         & m => this % mass, &
         & r1m => this % first_moment_mass, &
         & r2m => this % second_moment_mass )
      
      !residual = residual + m*qddot - m*g
      
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

    class(rigidbody)   , intent(inout) :: this
    type(scalar) , intent(inout) :: jacobian(:,:)
    type(scalar) , intent(in)    :: coeff(:)
    type(scalar) , intent(in)    :: U(:,:)
    
    associate( q => U(1,:), qdot => U(2,:), qddot => U(3,:), &
         & m => this % mass, &
         & r1m => this % first_moment_mass, &
         & r2m => this % second_moment_mass )
      
      !jacobian = jacobian + gamma*M
      
    end associate

  end subroutine add_jacobian
  
  !===================================================================!
  ! Sets the initial condition for use in the integator. If first order
  ! system just set initial Q, if a second order system set initial Q
  ! and qdot
  !===================================================================!  

  pure subroutine get_initial_condition(this, U)
    
    class(rigidbody) , intent(in)    :: this
    type(scalar)     , intent(inout) :: U(:,:)
    
    U(1,1:3) = this % position ! location (might use displacement)
    U(1,4:6) = this % orientation ! orientation (might use ang. displacement)

    U(2,1:3) = this % linear_velocity ! velocity
    U(2,4:6) = this % angular_velocity ! angular velocity
        
  end subroutine get_initial_condition

end module rigidbody_dynamics_class
