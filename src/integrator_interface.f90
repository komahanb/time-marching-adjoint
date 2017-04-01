#include "scalar.fpp"

!=====================================================================!
! Parent class for integration schemes to extend.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module integrator_interface

  use dynamic_physics_interface, only : dynamics

  implicit none

  private
  public ::  integrator

  type, abstract :: integrator

     !----------------------------------------------------------------!
     ! Contains the actual physical system
     !----------------------------------------------------------------!

     class(dynamics), pointer :: system => null()

     type(scalar)  :: tinit
     type(scalar)  :: tfinal
     type(scalar)  :: h

     !----------------------------------------------------------------!
     ! Track global time and states
     !----------------------------------------------------------------!

     type(scalar), allocatable :: time (:) ! time values (steps)
     type(scalar), allocatable :: U(:,:,:) ! state varibles (steps, deriv_ord, nvars)

     !----------------------------------------------------------------!
     ! Variables for managing time marching
     !----------------------------------------------------------------!

     type(logical) :: implicit
     type(integer) :: num_stages
     type(integer) :: num_steps
     type(integer) :: num_variables
     type(integer) :: time_deriv_order
     type(integer) :: print_level
     type(logical) :: approximate_jacobian

   contains

     procedure :: construct, destruct
     
     procedure :: evaluate_time

     !----------------------------------------------------------------!
     ! Deferred procedures for subtypes to implement                  !
     !----------------------------------------------------------------!

     procedure(evaluate_states_interface), deferred :: evaluate_states

     !----------------------------------------------------------------!
     ! Procedures                                                     !
     !----------------------------------------------------------------!

     procedure :: get_num_stages    , set_num_stages
     procedure :: get_num_steps     , set_num_steps
     procedure :: get_num_variables , set_num_variables
     procedure :: is_implicit       , set_implicit
     procedure :: set_physics
     procedure :: set_print_level
     procedure :: set_approximate_jacobian
     
     procedure :: integrate
     procedure :: write_solution
     procedure :: to_string
     
  end type integrator

  ! Define interfaces to deferred procedures
  interface

     !================================================================!
     ! Interface to approximate states using the time marching coeffs
     !================================================================!
     
     impure subroutine evaluate_states_interface(this, unew, uold)

       import integrator

       class(integrator), intent(in)  :: this
       type(scalar)     , intent(in)  :: uold(:,:,:)  ! previous values of state variables
       type(scalar)     , intent(out) :: unew(:,:)    ! approximated value at current step

     end subroutine evaluate_states_interface
     
  end interface

contains
   
  !===================================================================!
  ! Evaluate the indepedent variable (time)
  !===================================================================!
  
  impure subroutine evaluate_time(this, tnew, told, h)

    class(integrator) , intent(in)  :: this
    type(scalar)      , intent(in)  :: told    ! previous value of time
    type(scalar)      , intent(in)  :: h       ! step size
    type(scalar)      , intent(out) :: tnew    ! current time value
    
    advance_time: block
      
      type(scalar) :: tcoeff
      
      tnew   = told +  h
      
    end block advance_time

  end subroutine evaluate_time
    
  !===================================================================!
  ! Base class constructor logic
  !===================================================================!

  subroutine construct(this, system, tinit, tfinal, h, implicit)

    class(integrator) , intent(inout)         :: this
    class(dynamics)   , intent(in)   , target :: system
    type(scalar)      , intent(in)            :: tinit, tfinal
    type(scalar)      , intent(in)            :: h
    type(logical)     , intent(in)            :: implicit
    
    call this % set_physics(system)    

    this % tinit = tinit
    this % tfinal = tfinal
    this % h = h 

    call this % set_implicit(implicit)

    this % time_deriv_order = system % get_time_deriv_order()

    this % num_steps = int((this % tfinal - this % tinit)/this % h) + 1

    call this % set_num_variables( this % system % get_num_state_vars() )

    if ( .not. (this % get_num_variables() .gt. 0) ) then
       stop ">> Error: No DOF to march in time. Stopping."
    end if

    allocate(this % time( this % get_num_steps() ))
    this % time = 0.0d0
       
    allocate( this % U( &
         & this % get_num_steps(), &
         & this % time_deriv_order + 1, &
         & this % get_num_variables() &
         & ))
    this % U = 0.0d0

  end subroutine construct

  !======================================================================!
  ! Base class destructor
  !======================================================================!

  impure subroutine destruct(this)

    class(integrator), intent(inout) :: this

    ! Clear global states and time
    if(allocated(this % U)) deallocate(this % U)
    if(allocated(this % time)) deallocate(this % time)

  end subroutine destruct

  !===================================================================!
  ! Write solution to file
  !===================================================================!

  subroutine write_solution(this, filename)

    class(integrator)             :: this
    character(len=*), intent(in)  :: filename
    character(len=7), parameter   :: directory = "output/"
    character(len=:), allocatable :: path
    character(len=:), allocatable :: new_name
    integer                       :: k, j, i, ierr
    integer                       :: nsteps

    ! Open resource
    path = trim(filename)

    open(unit=90, file=trim(path), iostat= ierr)
    if (ierr .ne. 0) then
       write(*,'("  >> Opening file ", 39A, " failed")') path
       return
    end if
    
    ! Write data
    nsteps = this % get_num_steps()
    loop_time: do k = 1, nsteps
       write(90, *)  this % time(k), this % U (k,:,:)
    end do loop_time
    
    ! Close resource
    close(90)

  end subroutine write_solution

  !===================================================================!
  ! Time integration logic
  !===================================================================!

  impure subroutine integrate( this )    
  
    class(integrator), intent(inout) :: this
    integer :: k

    ! Set states to zero
    this % U     = 0.0d0
    this % time  = 0.0d0

    ! Get the initial condition
    call this % system % get_initial_condition(this % U(1,:,:))
    
    ! March in time
    time: do k = 2, this % num_steps

       call this % evaluate_time(this % time(k), this % time(k-1), this % h)

       call this % evaluate_states(this % U(k,:,:), this % U(1:k-1,:,:))

    end do time

  end subroutine integrate
  
  !===================================================================!
  ! Returns the number of stages per time step
  !===================================================================!
  
  impure type(integer) function get_num_stages(this)

    class(integrator), intent(in) :: this

    get_num_stages = this % num_stages

  end function get_num_stages

  !===================================================================!
  ! Sets the number of stages per time step
  !===================================================================!

  impure subroutine set_num_stages(this, num_stages)

    class(integrator), intent(inout) :: this
    type(integer)    , intent(in)    :: num_stages

    this % num_stages = num_stages

  end subroutine set_num_stages

  !===================================================================!
  ! Returns the number of steps
  !===================================================================!

  impure type(integer) function get_num_steps(this)

    class(integrator), intent(in) :: this

    get_num_steps = this % num_steps

  end function get_num_steps

  !===================================================================!
  ! Sets the number of steps
  !===================================================================!

  impure subroutine set_num_steps(this, num_steps)

    class(integrator), intent(inout) :: this
    type(integer)    , intent(in)    :: num_steps

    this % num_steps = num_steps

  end subroutine set_num_steps

  !===================================================================!
  ! Returns the number of variables (dof)
  !===================================================================!

  impure type(integer) function get_num_variables(this)

    class(integrator), intent(in) :: this

    get_num_variables = this % num_variables

  end function get_num_variables

  !===================================================================!
  ! Sets the number of variables (dof)
  !===================================================================!

  impure subroutine set_num_variables(this, num_variables)

    class(integrator), intent(inout) :: this
    type(integer)    , intent(in)    :: num_variables

    this % num_variables = num_variables

  end subroutine set_num_variables

  !===================================================================!
  ! See if the scheme is implicit
  !===================================================================!

  impure type(logical) function is_implicit(this)

    class(integrator), intent(in) :: this

    is_implicit = this % implicit

  end function is_implicit
  
  !===================================================================!
  ! Sets the scheme as implicit
  !===================================================================!
  
  impure subroutine set_implicit(this, implicit)

    class(integrator), intent(inout) :: this
    type(logical)    , intent(in)    :: implicit

    this % implicit = implicit

  end subroutine set_implicit

  !===================================================================!
  ! Setter that can be used to set the method in which jacobian needs
  ! to be computed. Setting this to .true. would make the code use
  ! finite differences, this is enabled by default too. If set to
  ! .false. the expects to provide implementation in assembleJacobian
  ! in a type that extends PHYSICS.
  !===================================================================!

  impure subroutine set_approximate_jacobian(this, approximate_jacobian)

    class(integrator), intent(inout) :: this
    type(logical), intent(in)     :: approximate_jacobian

    this % approximate_jacobian = approximate_jacobian

  end subroutine set_approximate_jacobian

  !===================================================================!
  ! Set ANY physical system that extends the type PHYSICS and provides
  ! implementation to the mandatory functions assembleResidual and
  ! getInitialStates
  !===================================================================!
  
  subroutine set_physics(this, physical_system)

    class(integrator)   :: this
    class(dynamics), target :: physical_system

    this % system => physical_system

  end subroutine set_physics

  !===================================================================!
  ! Manages the amount of print
  !===================================================================!

  impure subroutine set_print_level(this,print_level)

    class(integrator), intent(inout) :: this
    type(integer), intent(in)    :: print_level

    this % print_level = print_level

  end subroutine set_print_level

  !===================================================================!
  ! Prints important fields of the class
  !===================================================================!
  
  subroutine to_string(this)
    
    class(integrator), intent(in) :: this
    
    print '("  >> Physical System     : " ,A10)' , this % system % get_description()
    print '("  >> Start time          : " ,F8.3)', this % tinit
    print '("  >> End time            : " ,F8.3)', this % tfinal
    print '("  >> Step size           : " ,E9.3)', this % h
    print '("  >> Number of variables : " ,i4)'  , this % get_num_variables()
    print '("  >> Equation order      : " ,i4)'  , this % time_deriv_order
    print '("  >> Number of steps     : " ,i10)' , this % num_steps

  end subroutine to_string
  
end module integrator_interface
