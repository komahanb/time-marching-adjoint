#include "scalar.fpp"
!=====================================================================!
! Module that provides the ability to the user to implement random
! ODEs and use with the DIRK scheme
!=====================================================================!

Module vanderpol_class

!  use iso_fortran_env , only : dp => real64

  use physics_class,  only : physics
  use function_class, only : abstract_function
  
  implicit none

  private

  public :: vanderpol

  !-------------------------------------------------------------------!
  ! Type that implements vanderpol equations in first order form
  !-------------------------------------------------------------------!
  
  type, extends(physics) :: vanderpol

     ! Define constants and other parameters needed for residual and
     ! jacobian assembly here

     type(scalar) :: m = 1.0d0

   contains
     procedure :: mapDesignVars
     procedure :: assembleResidual
     procedure :: assembleJacobian
     procedure :: getInitialStates
     procedure :: getResidualDVSens

  end type vanderpol

contains
  
  !-------------------------------------------------------------------!
  ! Map the the design variables into the class variables
  !-------------------------------------------------------------------!
  
  subroutine mapDesignVars(this)

    class(vanderpol) :: this

    this % m = this % x(1)

  end subroutine mapDesignVars

!!$  !===================================================================!
!!$  ! Sets the design variables into the system
!!$  !===================================================================!
!!$  
!!$  subroutine setDesignVars(this, x)
!!$
!!$    class(vanderpol)                   :: this
!!$    real(8), intent(in), dimension(:)  :: x
!!$
!!$    ! Overwrite the values to supplied ones
!!$    if (this % num_design_vars .eq. 1) then 
!!$
!!$    else if (this % num_design_vars .eq. 2) then
!!$
!!$    else if (this % num_design_vars .eq. 3) then
!!$
!!$    end if
!!$
!!$  end subroutine setDesignVars

  !-------------------------------------------------------------------!
  ! Residual assembly at each time step. This is a mandary function
  ! that the user needs to implement to use the integration
  ! scheme. This is where the differential equations are supplied to
  ! the solver.
  ! -------------------------------------------------------------------!
  
  subroutine assembleResidual( this, res, time, u, udot, uddot )

    class(vanderpol) :: this
    type(scalar), intent(inout), dimension(:) :: res
    real(dp), intent(in)                      :: time
    type(scalar), intent(in), dimension(:)    :: u, udot, uddot

    res(1) = udot(1) - u(2)
    res(2) = udot(2) - this % m *( 1.0d0 - u(1)*u(1) )*u(2) + u(1)

  end subroutine assembleResidual

  !-------------------------------------------------------------------!
  ! Jacobian assembly at each time step. If you don't provide the
  ! analytical jacobian, set setApproximateJacobian(.true.) into the
  ! integrator object. We use finite-difference method to approximate
  ! the Jacobian.
  !
  ! Jacobian is the matrix of partial derivatives. Each row in the
  ! Jacobian matrix arises from differntiating a single equation. Each
  ! column in the Jacobian comes from a variable in the problem. Note
  ! the NEQN should be equal to NVARS for the system to be solved.
  !
  ! Note: alpha, beta and gamma are scalars that need to be multiplied
  ! with the partial derivatives DRDQ, DRDudot and DRDuddot
  ! respectively.
  ! -------------------------------------------------------------------!

  subroutine assembleJacobian( this, jac, alpha, beta, gamma, &
       & time, u, udot, uddot )

    class(vanderpol) :: this
    type(scalar), intent(inout), dimension(:,:) :: jac
    type(scalar), intent(in)                    :: alpha, beta, gamma
    real(dp), intent(in)                        :: time
    type(scalar), intent(in), dimension(:)      :: u, udot, uddot

    ! Zero all entries first
    jac = 0.0d0

    !-----------------------------------------------------------------!
    ! Add dR/dQ
    !-----------------------------------------------------------------!

    ! derivative of first equation
    
    jac(1,1) = jac(1,1) + alpha*0.0d0
    jac(1,2) = jac(1,2) - alpha*1.0d0

    ! derivative of second equation
    
    jac(2,1) = jac(2,1) + this % m*alpha*(1.0d0 + 2.0d0*u(1)*u(2))
    jac(2,2) = jac(2,2) + this % m*alpha*(u(1)*u(1)-1.0d0)
    
    !-----------------------------------------------------------------!
    ! Add dR/dQDOT
    !-----------------------------------------------------------------!
    
    ! derivative of first equation
    
    jac(1,1) = jac(1,1) + beta*1.0d0
    jac(1,2) = jac(1,2) + beta*0.0d0

    ! derivative of second equation

    jac(2,1) = jac(2,1) + this % m*beta*0.0d0
    jac(2,2) = jac(2,2) + this % m*beta*1.0d0

  end subroutine assembleJacobian

  !---------------------------------------------------------------------!
  ! Sets the initial condition for use in the integator. If first order
  ! system just set initial Q, if a second order system set initial Q
  ! and udot
  !---------------------------------------------------------------------!

  subroutine getInitialStates(this, time, u, udot)

    class(vanderpol) :: this


    real(dp), intent(in) :: time
    type(scalar), intent(inout), dimension(:) :: u, udot
    
    u(1) = 2.0d0
    u(2) = 0.0d0

  end subroutine getInitialStates

!!$  !===================================================================!
!!$  ! Return the number of state variables
!!$  !===================================================================!
!!$  
!!$  function getNumStateVars(this)
!!$
!!$    class(vanderpol) :: this
!!$    integer          :: getNumStateVars
!!$
!!$    getNumStateVars = this % num_state_vars
!!$
!!$  end function getNumStateVars

  !-------------------------------------------------------------------!
  ! Routine for evaluating the gradient of Residual with respect
  ! to the design X
  !-------------------------------------------------------------------!
  
  subroutine getResidualDVSens(this, jac, scale, time, x, u, udot, uddot)

    class(vanderpol)                       :: this
    type(scalar), intent(inout), dimension(:,:) :: jac
    real(dp), intent(in) :: time
    type(scalar), intent(in), dimension(:)      :: x, u, udot, uddot
    type(scalar)                                :: scale

    stop"Not implemented"

  end subroutine getResidualDVSens

end module vanderpol_class

