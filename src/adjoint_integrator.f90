#include "scalar.fpp"

module adjoint_integrator_interface

  use integrator_interface, only : integrator

  implicit none
  
  private
  public :: adjoint_integrator
  
  type, abstract, extends(integrator) :: adjoint_integrator

     ! Takes an instance of integrator to solve for state variables
     class(integrator), allocatable :: iforward
          
     ! Need to have a function

   contains

     ! Overrider the solution procedure tailored for adjoint solve
     procedure :: solve => adjoint_solve

  end type adjoint_integrator

contains

  !==================================================================!
  ! Solving the adjoint linear system
  !==================================================================!
  
  subroutine find_adjoint(system, coeff, t, Q)

    use dynamic_physics_interface, only : dynamics
    use nonlinear_algebra, only : approximate_jacobian
    use linear_algebra

    class(dynamics) , intent(inout) :: system
    type(scalar)    , intent(in)    :: coeff(:)
    type(scalar)    , intent(in)    :: t
    type(scalar)    , intent(inout) :: Q(:,:)

    ! Other Local variables
    type(scalar), allocatable, dimension(:)   :: res
    type(scalar), allocatable, dimension(:,:) :: jac, fd_jac
    integer                                   :: nvars
    type(scalar)                              :: jac_err
    type(scalar), allocatable                 :: U(:,:)

    ! 
    type(logical) :: jacobian_check = .true.

    ! find the size of the linear system based on the calling object
    nvars = size(Q(1,:))

    if ( nvars .ne. system % get_num_state_vars() ) stop "NVARS-MISMATCH"

    if ( .not. allocated(res)    ) allocate( res(nvars)          )
    if ( .not. allocated(jac)    ) allocate( jac(nvars,nvars)    )
    if ( .not. allocated(fd_jac) ) allocate( fd_jac(nvars,nvars) )

    res = 0.0d0
    call system % add_residual(res, Q)

    ! Get the jacobian matrix
    if ( system % is_approximate_jacobian() ) then

       ! Compute an approximate Jacobian using finite differences
       allocate(U, source=Q)
       jac = 0.0d0
       call approximate_jacobian(system, jac, coeff, U)
       deallocate(U)

    else

       ! Use the user supplied Jacobian implementation
       jac = 0.0d0
       call system % add_jacobian(jac, coeff, Q)

       ! Check the Jacobian implementation once at the beginning of integration
       if (   system % get_differential_order() .gt. 0 .and. &
            & system % get_differential_order() .lt. 3 .and. &
            & jacobian_check ) then

          ! Compute an approximate Jacobian using finite differences
          allocate(U, source=Q)
          call approximate_jacobian(system, fd_jac, coeff, U)
          deallocate(U)

          ! Compare the exact and approximate Jacobians and
          ! complain about the error in Jacobian if there is any
          jac_err = maxval(abs(fd_jac - jac))
          if ( abs(jac_err) .gt. 1.0d-6) then
             print *, "q     =", Q(1,:)
             print *, "qdot  =", Q(2,:)
             print *, "qddot =", Q(3,:)
             print *, "a,b,c =", coeff
             print *, "J     =", jac
             print *, "Jhat  =", fd_jac
             print *, "WARNING: Possible error in function jacobian", jac_err
          end if

          ! Set that the jacobian is checked
          jacobian_check = .false.

       end if

    end if

    ! Call LAPACK to solve the linear system with transposed jacobian
    jac = transpose(jac)
    U(1,:) = solve(jac, -res)

    if (allocated(res))    deallocate(res)
    if (allocated(jac))    deallocate(jac)
    if (allocated(fd_jac)) deallocate(fd_jac)

  end subroutine find_adjoint
  
  !===================================================================!
  ! Time integration logic
  !===================================================================!
  
  impure subroutine adjoint_solve( this )

    class(adjoint_integrator), intent(inout) :: this

    ! Perform a forward solution to find the state variables in time
    call this % iforward % solve()

    ! Perform a reverse solve to find the adjoint variables in time
    
    
  end subroutine adjoint_solve
  
end module adjoint_integrator_interface

!=====================================================================!
! Implementation of BDF discrete adjoint
!=====================================================================!

module adjoint_backward_differences_class

  use integrator_interface, only : integrator
  use adjoint_integrator_interface, only : adjoint_integrator

  implicit none
  
  private
  public :: BDFDA
  
  ! Constructor
  interface BDFDA
     module procedure create
  end interface BDFDA
  
  ! Type definition for BDF discrete adjoint
  type, extends(adjoint_integrator) :: BDFDA

   contains

     ! Implement deferred procedures from base class
     procedure :: step => adjoint_step
     procedure :: get_bandwidth => adjoint_bandwidth
     
     ! Destructor
     final :: destroy

  end type BDFDA

contains

  !===================================================================!
  ! Initialize the BDF Adjoint datatype and allocate required
  ! variables
  !===================================================================!
  
  type(bdfda) function create(forward) result(this)

    class(integrator), intent(in) :: forward

    ! Store the forward integration object 
    allocate(this % iforward, source = forward)    

    ! Set the rest of the attributes
    call this % construct(      &
         & forward % system,    &
         & forward % tinit,     &
         & forward % tfinal,    &
         & forward % h,         &
         & forward % implicit,  &
         & forward % num_stages )

  end function create

  !===================================================================!
  ! Destructor
  !===================================================================!
  
  pure subroutine destroy(this)

    type(BDFDA), intent(inout) :: this

    ! Parent class call
    call this % destruct()

    ! Deallocate the forward integrator
    if(allocated(this % iforward)) deallocate(this % iforward)

  end subroutine destroy

  !===================================================================!
  ! Returns the order of approximation for the given time step k and
  ! degree d
  !===================================================================!

  pure type(integer) function adjoint_bandwidth(this, time_index) result(width)
    
    class(BDFDA), intent(in) :: this
    type(integer), intent(in) :: time_index

    width = this % iforward % get_bandwidth(time_index)
    
  end function adjoint_bandwidth

  !===================================================================!
  ! Take a time step using the supplied time step and order of
  ! accuracy
  !===================================================================!

  impure subroutine adjoint_step(this, t, u, h, p, ierr)

    use linear_algebra, only : solve

    ! Argument variables
    class(bdfda)  , intent(inout) :: this
    type(scalar)  , intent(inout) :: t(:)
    type(scalar)  , intent(inout) :: u(:,:,:)
    type(integer) , intent(in)    :: p
    type(scalar)  , intent(in)    :: h
    type(integer) , intent(out)   :: ierr

    ! Local variables
    type(scalar) , allocatable :: lincoeff(:) ! order of equation + 1
    type(integer) :: torder, n, i, k
    type(scalar)  :: scale

    associate(psi => U)
      
    end associate

  end subroutine adjoint_step

end module adjoint_backward_differences_class
