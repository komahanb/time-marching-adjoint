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
