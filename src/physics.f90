#include "scalar.fpp"
!=====================================================================!
! Module that contains common procedures for any physical system
! subject to governing equations
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module physics_class

  use function_class, only  : abstract_function

  implicit none

  private
  public :: physics
  
  !-------------------------------------------------------------------!
  ! Type that models any physical phenomenon
  !-------------------------------------------------------------------!

  type, abstract :: physics

     integer                             :: num_design_vars = 0 
     integer                             :: num_state_vars  = 0 
     type(scalar), dimension(:), allocatable :: x
     class(abstract_function), pointer   :: func => null() ! function of interest
     
   contains  

     ! Constructor and destructor
     procedure :: initialize, finalize

     ! Setters
     procedure :: setFunction
     procedure :: setDesignVars
     procedure :: getNumStateVars

     ! Deferred procedures
     procedure(initial_condition_interface), deferred :: getInitialStates     
     procedure(residual_assembly_interface), deferred :: assembleResidual
     procedure(jacobian_assembly_interface), deferred :: assembleJacobian
 
     procedure(InterfaceGetResidualDVSens), deferred  :: getResidualDVSens
     procedure(InterfaceMapDesignVars), deferred  :: mapDesignVars

  end type physics
  
  interface
     
     !----------------------------------------------------------------!
     ! Interface for evaluating the gradient of Residual with respect
     ! to X
     !----------------------------------------------------------------!
     
     subroutine InterfaceGetResidualDVSens(this, jac, scale, time, x, u, udot, uddot)

       import physics
       
       class(physics)                         :: this
       type(scalar), intent(inout), dimension(:,:) :: jac
       real(dp), intent(in)                        :: time
       type(scalar), intent(in), dimension(:)      :: x, u, udot, uddot
       type(scalar)                                :: scale

     end subroutine InterfaceGetResidualDVSens

     !----------------------------------------------------------------!
     ! Interface for initialization tasks
     !----------------------------------------------------------------!
     
     subroutine InterfaceInitialize(this,  x, function)
       import physics
       import abstract_function
       class(physics) :: this
       class(abstract_function), target, OPTIONAL  :: function
       type(scalar), intent(in), dimension(:), OPTIONAl :: x
     end subroutine InterfaceInitialize
    
     !----------------------------------------------------------------!
     ! Interface for residual assembly at each time step
     !----------------------------------------------------------------!

     subroutine residual_assembly_interface(this, res, time, u, udot, uddot)

       import physics

       class(physics) :: this
       type(scalar), intent(inout), dimension(:) :: res
       real(dp), intent(in) :: time
       type(scalar), intent(in), dimension(:) :: u, udot, uddot

     end subroutine residual_assembly_interface

     !----------------------------------------------------------------!
     ! Interface for jacobian assembly at each time step
     !----------------------------------------------------------------!

     subroutine jacobian_assembly_interface(this, jac, alpha, beta, gamma, &
          & time, u, udot, uddot)

       import physics

       class(physics) :: this
       type(scalar), intent(inout), dimension(:,:) :: jac
       type(scalar), intent(in) :: alpha, beta, gamma
       real(dp), intent(in) :: time
       type(scalar), intent(in), dimension(:) :: u, udot, uddot

     end subroutine jacobian_assembly_interface
     
     !----------------------------------------------------------------!
     ! Interface for supplying the initial condition to the integrator!
     !----------------------------------------------------------------!
     
     subroutine initial_condition_interface(this, time, u, udot)

       import physics

       class(physics) :: this
       real(dp), intent(in) :: time
       type(scalar), intent(inout), dimension(:) :: u, udot

     end subroutine initial_condition_interface

     !----------------------------------------------------------------!
     ! Return the number of state variables
     !----------------------------------------------------------------!
     
     function InterfaceGetNumStateVars(this)
       import physics
       class(physics) :: this
       integer :: InterfaceGetNumStateVars
     end function InterfaceGetNumStateVars
     
     !-------------------------------------------------------------------!
     ! Map the the design variables into the class variables
     !-------------------------------------------------------------------!
     
     subroutine InterfaceMapDesignVars(this)
       import physics
       class(physics) :: this
     end subroutine InterfaceMapDesignVars

  end interface

contains

  !===================================================================!
  ! Initialize the system
  !===================================================================!

  subroutine initialize(this, num_state_vars, num_design_vars)

    class(physics)      :: this
    integer, intent(in) :: num_state_vars
    integer, intent(in), OPTIONAL :: num_design_vars
    
    this % num_state_vars  = num_state_vars
    
    if (present(num_design_vars)) then
       this % num_design_vars = num_design_vars
       allocate(this % x(this % num_design_vars))
    end if

  end subroutine initialize

  !===================================================================!
  ! Finalize the allocated variables
  !===================================================================!

  subroutine finalize(this)

    class(physics) :: this
    
    if ( allocated(this % x) ) deallocate(this % x)

    print*, "Yet to deassociate the function"

  end subroutine finalize

  !-------------------------------------------------------------------!
  ! Set the design varaibles
  !-------------------------------------------------------------------!
  
  subroutine setDesignVars(this, num_dvs, x)
    
    class(physics)                    :: this
    type(scalar), intent(in), dimension(:) :: x
    integer, intent(in)               :: num_dvs
    
    if (this % num_design_vars .ne. size(x)) stop "Error in num_design_vars"

    this % x = x

    call this % mapDesignVars()

  end subroutine setDesignVars
    
  !===================================================================!
  ! Set the function created into the system                          !
  !===================================================================!
  
  subroutine setFunction(this, func)

    class(physics)                   :: this
    class(abstract_function), target :: func
    
    this % func => func
    
  end subroutine setFunction
  
  !===================================================================!
  ! Return the number of state variables
  !===================================================================!
  
  integer function getNumStateVars(this)

    class(physics) :: this

    getNumStateVars = this % num_state_vars

  end function getNumStateVars

end module physics_class
