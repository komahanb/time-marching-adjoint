#include "scalar.fpp"

!=====================================================================!
! VANDERPOL OSCILLATOR physics. This is a simple example of setting up
! a physical system for analysis within the framework.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module vanderpol_system

  use constants                        , only : WP
  use time_dependent_physics_interface , only : dynamics

  ! Currently will use dense matrix implementaion 
  use dense_vector_interface           , only : dense_vector
  use dense_matrix_interface           , only : dense_matrix

  use vector_interface           , only : vector
  use matrix_interface           , only : matrix

  implicit none

  private

  public :: vanderpol

  !-------------------------------------------------------------------!
  ! Type that implements vanderpol equations in first order form
  !-------------------------------------------------------------------!
  
  type, extends(dynamics) :: vanderpol

     ! Define constants and other parameters needed for residual and
     ! jacobian assembly here

     type(scalar) :: m = 1.0_WP

   contains

     ! Implement deferred procedures from superclass
     procedure :: assemble_residual
     procedure :: assemble_jacobian
     procedure :: get_initial_condition

  end type vanderpol

  ! Interface to construct vanderpol system
  interface vanderpol
     procedure constructor
  end interface vanderpol

contains
 
  !===================================================================!
  ! Constructor for dense matrix
  !===================================================================!
  
  pure type(vanderpol) function constructor(mass) result (this)

    type(scalar), intent(in) :: mass

    ! Set time order of vanderpol system
    call this % set_time_order(1)

    ! Set the oscillator parameter
    this % m = mass

  end function constructor

  !-------------------------------------------------------------------!
  ! Residual assembly at each time step. This is a mandary function
  ! that the user needs to implement to use the integration
  ! scheme. This is where the differential equations are supplied to
  ! the solver.
  ! -------------------------------------------------------------------!
  
  pure subroutine assemble_residual(this, residual, state_vectors)

    class(vanderpol), intent(inout) :: this
    class(vector)   , intent(inout) :: residual
    class(vector)   , intent(in)    :: state_vectors(:)

!!$    select type(residual)
!!$    class is (dense_vector)
!!$       call assemble_dense_residual( residual % vals )
!!$    class default
!!$    end select

  end subroutine assemble_residual
  
!!$  pure subroutine assemble_dense_residual(R, )
!!$
!!$    set_vales: block
!!$      type(scalar) :: r( residual % get_size() )
!!$      call residual % set_entry(1, 
!!$      r(1) = udot(1) - u(2)
!!$      r(2) = udot(2) - this % m*(1.0_wp-u(1)*u(1))*u(2) + u(1)
!!$    end block set_vales
!!$
!!$  end subroutine assemble_dense_residual

  !-------------------------------------------------------------------!
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
  ! -------------------------------------------------------------------!
  
  pure subroutine assemble_jacobian(this, jacobian, state_vectors, coeffs)
    
    class(vanderpol) , intent(inout) :: this
    class(matrix)    , intent(inout) :: jacobian
    class(vector)    , intent(in)    :: state_vectors(:)
    type(scalar)     , intent(in)    :: coeffs(:)
    
!!$    associate( u => state_vectors(1) % vals, &
!!$         & udot => state_vectors(2) % vals, &
!!$         & jac => jacobian % vals, &
!!$         & alpha => coeffs(1), &
!!$         & beta => coeffs(2) )
!!$
!!$      ! Zero all entries first
!!$      jac = 0.0_WP
!!$
!!$      !-----------------------------------------------------------------!
!!$      ! Add dR/dQ
!!$      !-----------------------------------------------------------------!
!!$
!!$      ! derivative of first equation
!!$
!!$      jac(1,1) = jac(1,1) + alpha*0.0_WP
!!$      jac(1,2) = jac(1,2) - alpha*1.0_WP
!!$
!!$      ! derivative of second equation
!!$
!!$      jac(2,1) = jac(2,1) + this % m*alpha*(1.0_WP + 2.0_WP*u(1)*u(2))
!!$      jac(2,2) = jac(2,2) + this % m*alpha*(u(1)*u(1)-1.0_WP)
!!$
!!$      !-----------------------------------------------------------------!
!!$      ! Add dR/dQDOT
!!$      !-----------------------------------------------------------------!
!!$
!!$      ! derivative of first equation
!!$
!!$      jac(1,1) = jac(1,1) + beta*1.0_WP
!!$      jac(1,2) = jac(1,2) + beta*0.0_WP
!!$
!!$      ! derivative of second equation
!!$
!!$      jac(2,1) = jac(2,1) + this % m*beta*0.0_WP
!!$      jac(2,2) = jac(2,2) + this % m*beta*1.0_WP
!!$
!!$    end associate

  end subroutine assemble_jacobian
  
  !---------------------------------------------------------------------!
  ! Sets the initial condition for use in the integator. If first order
  ! system just set initial Q, if a second order system set initial Q
  ! and udot
  !---------------------------------------------------------------------!

  pure subroutine get_initial_condition(this, state_vectors)

    class(vanderpol) , intent(inout) :: this
    class(vector)    , intent(inout) :: state_vectors(:)

!!$    associate( u => state_vectors(1) % vals, &
!!$         & udot => state_vectors(2) % vals )
!!$
!!$      u(1) = 2.0_WP
!!$      u(2) = 0.0_WP
!!$
!!$    end associate

  end subroutine get_initial_condition

end module vanderpol_system
