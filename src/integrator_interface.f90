#include "scalar.fpp"

!=====================================================================!
! Parent class for integration schemes to extend. This has some common
! logic such as:
!
! (1) nonlinear solution process
! (2) approximating derivatives 
! (3) writing solution to files
! (4) adjoint system solution
!
!=====================================================================!

module integrator_class

  use physics_interface, only : physics

  implicit none

  private
  public ::  integrator

  !===================================================================! 
  ! Abstract Integrator type
  !===================================================================! 

  type, abstract :: integrator

     !----------------------------------------------------------------!
     ! Contains the actual physical system
     !----------------------------------------------------------------!

     class(physics), pointer :: system => null()

     type(scalar)  :: tinit     = 0.0d0 ! initial time
     type(scalar)  :: tfinal    = 1.0d0 ! final time
     type(scalar)  :: h         = 0.1d0 ! default step size

     !----------------------------------------------------------------!
     ! Track global time and states
     !----------------------------------------------------------------!

     type(scalar), allocatable :: time (:)     ! time values (steps)
     type(scalar), allocatable :: U    (:,:,:) ! state varibles (steps, eqn_ord, nvars)

     !----------------------------------------------------------------!
     ! Variables for managing time marching
     !----------------------------------------------------------------!

     type(logical) :: implicit             = .true.
     type(integer) :: num_stages           = 0
     type(integer) :: num_steps            = 0
     type(integer) :: num_variables        = 0
     type(integer) :: equation_order       = 0
     type(integer) :: print_level          = 0
     type(logical) :: approximate_jacobian = .false.

   contains

     !----------------------------------------------------------------!
     ! Deferred procedures for subtypes to implement                  !
     !----------------------------------------------------------------!

     procedure(evaluate_states_interface)        , private, deferred :: evaluate_states
     procedure(get_state_coefficients_interface) , private, deferred :: get_state_coefficients
     procedure(get_state_coefficients_interface) , private, deferred :: get_time_coefficients

     !----------------------------------------------------------------!
     ! Procedures                                                     !
     !----------------------------------------------------------------!

     procedure :: get_num_stages    , set_num_stages
     procedure :: get_num_steps     , set_num_steps
     procedure :: get_num_variables , set_num_variables
     procedure :: is_implicit       , set_implicit
     
     procedure :: evaluate_time
     procedure :: writeSolution
     procedure :: setPhysicalSystem 
     procedure :: setPrintLevel
     procedure :: setApproximateJacobian
     procedure :: construct, destruct
     procedure :: integrate
     
  end type integrator

  ! Define interfaces to deferred procedures
  abstract interface

     !================================================================!
     ! Interface to approximate states using the time marching coeffs
     !================================================================!
     
     subroutine evaluate_states_interface(this, unew, uold, scoeff, order, h)

       import integrator

       class(integrator), intent(inout) :: this
       type(scalar)     , intent(in)    :: uold(:,:)    ! previous values of state variables
       type(scalar)     , intent(out)   :: unew(:,:)    ! approximated value at current step
       type(scalar)     , intent(in)    :: scoeff(:,:)  ! alpha=scoeff(1,:), beta=scoeff(2,:), gamma=scoeff(3,:)
       type(integer)    , intent(in)    :: bound(:)     ! p
       type(scalar)     , intent(in)    :: h

       ! evaluate the state approximations based on the integration method

     end subroutine evaluate_states_interface

     !================================================================!
     ! Retrieve the state approximation coefficients
     !================================================================!
     
     subroutine get_state_coefficients_interface(this, coeff, int_order, eqn_order)

       import integrator

       class(integrator) , intent(inout) :: this
       type(integer)     , intent(in)    :: int_ord(:)  ! order of approximation of the integration
       type(integer)     , intent(in)    :: eqn_ord     ! order of the differential equation in time
       type(scalar)      , intent(out)   :: coeff(:)    ! order of equation + 1

       ! Set the coefficients in to the coeff array and return

     end subroutine get_state_coefficients_interface

  end interface

contains
  
  !===================================================================!
  ! Returns the number of stages per time step
  !===================================================================!
  
  pure type(integer) function get_num_stages(this)

    class(integrator), intent(in) :: this

    get_num_stages = this % num_stages

  end function get_num_stages

  !===================================================================!
  ! Sets the number of stages per time step
  !===================================================================!

  pure subroutine set_num_stages(this, num_stages)

    class(integrator), intent(inout) :: this
    type(integer)    , intent(in)    :: num_stages

    this % num_stages = num_stages

  end subroutine set_num_stages

  !===================================================================!
  ! Returns the number of steps
  !===================================================================!

  pure type(integer) function get_num_steps(this)

    class(integrator), intent(in) :: this

    get_num_steps = this % num_steps

  end function get_num_steps

  !===================================================================!
  ! Sets the number of steps
  !===================================================================!

  pure subroutine set_num_steps(this, num_steps)

    class(integrator), intent(inout) :: this
    type(integer)    , intent(in)    :: num_steps

    this % num_steps = num_steps

  end subroutine set_num_steps

  !===================================================================!
  ! Returns the number of variables (dof)
  !===================================================================!

  pure type(integer) function get_num_variables(this)

    class(integrator), intent(in) :: this

    get_num_variables = this % num_variables

  end function get_num_variables

  !===================================================================!
  ! Sets the number of variables (dof)
  !===================================================================!

  pure subroutine set_num_variables(this, num_variables)

    class(integrator), intent(inout) :: this
    type(integer)    , intent(in)    :: num_variables

    this % num_variables = num_variables

  end subroutine set_num_variables

  !===================================================================!
  ! See if the scheme is implicit
  !===================================================================!

  pure type(logical) function is_implicit(this)

    class(integrator), intent(in) :: this

    is_implicit = this % implicit

  end function is_implicit
  
  !===================================================================!
  ! Sets the scheme as implicit
  !===================================================================!
  
  pure subroutine set_implicit(this, implicit)

    class(integrator), intent(inout) :: this
    type(logical)    , intent(in)    :: implicit

    this % implicit = implicit

  end subroutine set_implicit

  !===================================================================!
  ! Evaluate the indepedent variable (time)
  !===================================================================!
  
  subroutine evaluate_time(this, tnew, told, tcoeff, h)

    class(integrator) , intent(inout) :: this
    type(scalar)      , intent(in)    :: told    ! previous value of time
    type(scalar)      , intent(in)    :: tcoeff  ! coeffcient to scale the timestep h
    type(scalar)      , intent(in)    :: h       ! step size
    type(scalar)      , intent(out)   :: tnew    ! current time value

    tnew = told + tcoeff * h

  end subroutine evaluate_time
  
  !===================================================================!
  ! Base class constructor logic
  !===================================================================!

  subroutine construct(this, system, tinit, tfinal, h)

    class(integrator) , intent(inout)          :: this
    class(physics)    , intent(in)    , target :: system
    real(dp)          , intent(in)             :: tinit, tfinal
    real(dp)          , intent(in)             :: h

    !-----------------------------------------------------------------!
    ! Set the physical system in to the integrator                    !
    !-----------------------------------------------------------------!

    call this % setPhysicalSystem(system)
    print '("  >> Physical System        : ",A10)', this % system % name

   !-----------------------------------------------------------------!
    ! Set the initial and final time
    !-----------------------------------------------------------------!

    if (present(tinit)) then
       this % tinit = tinit
    end if
    print '("  >> Start time             : ",F8.3)', this % tinit

    if (present(tfinal)) then
       this % tfinal = tfinal
    end if
    print '("  >> End time               : ",F8.3)', this % tfinal

    !-----------------------------------------------------------------!
    ! Set the user supplied initial step size
    !-----------------------------------------------------------------!

    if (present(h)) then
       this % h = h 
    end if
    print '("  >> Step size              : ",E9.3)', this % h

    !-----------------------------------------------------------------!
    ! Fetch the number of state variables from the system object
    !-----------------------------------------------------------------!

    call set_num_variables( system % get_num_vars )
    print '("  >> Number of variables    : ",i4)', this % get_num_variables()

    if ( .not. (this % get_num_variables() .gt. 0) ) then
       stop ">> Error: No state variable. Stopping."
    end if

    !-----------------------------------------------------------------!
    ! Set the order of the governing equations
    !-----------------------------------------------------------------!

    this % equation_order = system % get_time_order()
    print '("  >> Equation order           : ",i4)', this % equation_order

    !-----------------------------------------------------------------!
    ! Find the number of time steps required during integration
    !-----------------------------------------------------------------!

    this % num_steps = int((this % tfinal - this % tinit)/this % h) + 1 
    print '("  >> Number of steps        : ",i10)', this % num_steps

    !-----------------------------------------------------------------!
    ! Allocate space for the global states and time
    !-----------------------------------------------------------------!

    allocate(this % time( this % get_num_steps() ))
    this % time = 0.0d0
    this % time(1) = this % tinit

    allocate( this % U( &
         & this % get_num_steps(), &
         & this % equation_order, &
         & this % get_num_variables() &
         & ))
    this % U = 0.0d0

  end subroutine construct

  !======================================================================!
  ! Base class destructor
  !======================================================================!

  subroutine destruct(this)

    class(integrator) :: this

    ! Clear global states and time
    if(allocated(this % U)) deallocate(this % U)
    if(allocated(this % time)) deallocate(this % time)

  end subroutine destruct

  !===================================================================!
  ! Setter that can be used to set the method in which jacobian needs
  ! to be computed. Setting this to .true. would make the code use
  ! finite differences, this is enabled by default too. If set to
  ! .false. the expects to provide implementation in assembleJacobian
  ! in a type that extends PHYSICS.
  !===================================================================!

  subroutine setApproximateJacobian(this, approximateJacobian)

    class(integrator) :: this
    logical :: approximateJacobian

    this % approximate_jacobian = approximateJacobian

  end subroutine setApproximateJacobian

  !===================================================================!
  ! Set ANY physical system that extends the type PHYSICS and provides
  ! implementation to the mandatory functions assembleResidual and
  ! getInitialStates
  !===================================================================!

  subroutine setPhysicalSystem(this, physical_system)

    class(integrator)      :: this
    class(physics), target :: physical_system

    this % system => physical_system

  end subroutine setPhysicalSystem

  !===================================================================!
  ! Manages the amount of print
  !===================================================================!

  subroutine setPrintLevel(this,print_level)

    class(integrator) :: this
    integer           :: print_level

    this % print_level = print_level

  end subroutine setPrintLevel

  !===================================================================!
  ! Write solution to file
  !===================================================================!

  subroutine writeSolution(this, filename)

    class(integrator)                      :: this
    character(len=*), OPTIONAL, intent(in) :: filename
    character(len=7), parameter            :: directory = "output/"
    character(len=32)                      :: path = ""
    integer                                :: k, j, ierr

    path = trim(path)

    if (present(filename)) then
       path = directory//filename
    else
       path = directory//"solution.dat"
    end if

    open(unit=90, file=trim(path), iostat= ierr)

    if (ierr .ne. 0) then
       write(*,'("  >> Opening file ", 39A, " failed")') path
       return
    end if

    if (this % equation_order .eq. 2 ) then

       do k = 1, this % num_steps
          write(90, *)  this % time(k), &
               & (dble(this % U (k,1,j) ), j=1,this % get_num_variables()  ), &
               & (dble(this % U (k,2,j) ), j=1,this % get_num_variables()  ), &
               & (dble(this % U (k,3,j) ), j=1,this % get_num_variables()  )
       end do

    else

       do k = 1, this % num_steps
          write(90, *)  this % time(k), &
               & (dble(this % U (k,1,j) ), j=1,this % get_num_variables()  ), &
               & (dble(this % U (k,2,j) ), j=1,this % get_num_variables()  )
       end do

    end if

    close(90)

  end subroutine writeSolution

  !===================================================================!
  ! Time integration logic
  !===================================================================!

  subroutine integrate( this )
    
    use nonlinear_algebra, only: nonlinear_solve

    class(integrator) :: this
    type(scalar)      :: alpha, beta, gamma
    integer           :: k
    type(scalar)      :: coeff(3)

    ! Set states to zeror
    this % U     = 0.0d0
    this % UDOT  = 0.0d0
    this % UDDOT = 0.0d0
    this % time  = 0.0d0

    ! Set the initial condition
    call this % system % getInitialStates(this % time(1), &
         & this % u(1,:), this % udot(1,:))

    this % current_step = 1

    ! March in time
    time: do k = 2, this % num_steps

       this % current_step =  k

       ! Increment the time
       this % time(k) = this % time(k-1) + this % h

       ! Approximate the states u, udot and uddot using ABM stencil
       call this % approximateStates()

       ! Determine the coefficients for linearing the Residual
       call this % getLinearCoeff(alpha, beta, gamma)

       ! Solve the nonlinear system at each step by driving the
       ! residual to zero
       call nonlinear_solve(this % system, &
            & alpha, beta, gamma, &
            & this % time(k), this % u(k,:), this % udot(k,:), this % uddot(k,:))
       
    end do time

  end subroutine integrate
  
end module integrator_class
