#include "scalar.fpp"

!=====================================================================!
! Adams Bashworth Moulton Integration Module for n-th order
! differential systems.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================! 

module abm_integrator_class

  use integrator_interface      , only : integrator
  use dynamic_physics_interface , only : dynamics

  implicit none

  private
  public :: ABM
  
  !===================================================================! 
  ! ABM Integrator type
  !===================================================================! 

  type, extends(integrator) :: ABM
     
     private

     ! ABM variables
     type(integer) :: max_order = 6
     type(scalar), allocatable :: A(:,:)

   contains
           
     procedure :: get_linearization_coeff
     procedure :: step
     procedure :: get_bandwidth

     ! Destructor
     final :: destroy

  end type ABM

  interface ABM
     module procedure create
  end interface ABM

contains

  !===================================================================!
  ! Initialize the ABM datatype and allocate required variables
  !===================================================================!
  
  type(abm) function create(system, tinit, tfinal, h, implicit, accuracy_order) result(this)

    class(dynamics)   , intent(in)   , target :: system
    type(scalar)      , intent(in)            :: tinit, tfinal
    type(scalar)      , intent(in)            :: h
    type(integer)     , intent(in)            :: accuracy_order
    type(logical)     , intent(in)            :: implicit   

    print *, "======================================"
    print *, ">>   Adams Bashforth Moulton       << "
    print *, "======================================"

    call this % construct(system, tinit, tfinal, h, implicit, 0)

    if (accuracy_order .le. this % max_order) this % max_order = accuracy_order
    print '("  >> Max ABM Order          : ",i4)', this % max_order

    ! Set the coefficients
    allocate( this % A (this % max_order, this % max_order) )
    this % A = 0.0d0
    if ( this % max_order .ge. 1 ) this % A(1,1:1) = [1.0d0]
    if ( this % max_order .ge. 2 ) this % A(2,1:2) = [1.0d0, 1.0d0]/2.0d0
    if ( this % max_order .ge. 3 ) this % A(3,1:3) = [5.0d0, 8.0d0, -1.0d0]/12.0d0
    if ( this % max_order .ge. 4 ) this % A(4,1:4) = [9.0d0, 19.0d0, -5.0d0, 1.0d0]/24.0d0
    if ( this % max_order .ge. 5 ) this % A(5,1:5) = [251.0d0, 646.0d0, -264.0d0, 106.0d0, -19.0d0]/720.0d0
    if ( this % max_order .ge. 6 ) this % A(6,1:6) = [475.0d0, 1427.0d0, -798.0d0, 482.0d0, -173.0d0, 27.0d0]/1440.0d0

    ! Sanity check on ABM coeffs
    sanity_check: block
      type(integer) :: j
      do j = 1, this % max_order
         if ( real(sum(this % A(j,1:j)) - 1.0d0) .gt. 0.00001 ) then
            stop "Error in ABM Coeff"
         end if
      end do
    end block sanity_check

  end function create

  !=================================================================!
  ! Destructor for the ABM integrator
  !=================================================================!
  
  impure subroutine destroy(this)

    type(ABM), intent(inout) :: this

    ! Parent class call
    call this % destruct()

    ! Deallocate ABM coefficient
    if(allocated(this % A)) deallocate(this % A)

  end subroutine destroy

  !===================================================================!
  ! Returns the order of approximation for the given time step k and
  ! degree d
  !===================================================================!

  pure type(integer) function get_bandwidth(this, time_index) result(width)

    class(ABM)   , intent(in) :: this
    type(integer), intent(in) :: time_index

    width = time_index - 1

    if (width .gt. this % max_order) width = this % max_order

  end function get_bandwidth

  !================================================================!
  ! Take a time step using the supplied time step and order of
  ! accuracy
  ! ================================================================!

  impure subroutine step(this, t, u, h, p, ierr)

    use nonlinear_algebra, only : solve

    ! Argument variables
    class(ABM)   , intent(inout) :: this
    type(scalar) , intent(inout) :: t(:)
    type(scalar) , intent(inout) :: u(:,:,:)
    type(integer), intent(in)    :: p
    type(scalar) , intent(in)    :: h
    type(integer), intent(out)   :: ierr

    ! Local variables
    type(scalar), allocatable :: lincoeff(:)  ! order of equation + 1
    type(integer) :: torder, n, i, k
    type(scalar)  :: scale
    
    ! Determine remaining paramters needed
    k = size(u(:,1,1))
    torder = this % system % get_differential_order()

    ! Advance the time to next step
    t(k) = t(k-1) + h

    ! Assume a value for highest order state
    u(k,torder+1,:) = 0.0d0

    ! Find the lower order states based on ABM formula
    do n = torder, 1, -1
       u(k,n,:) = u(k-1,n,:)
       do i = 0, p-1
          scale = h*this % A(p,i+1) !(t(k-i)-t(k-i-1))
          u(k,n,:) = u(k,n,:) + scale*u(k-i,n+1,:)
       end do
    end do

    ! Perform a nonlinear solution if this is a implicit method
    if ( this % is_implicit() ) then
       allocate(lincoeff(torder+1))         
       call this % get_linearization_coeff(p, h, lincoeff)
       call solve(this % system, lincoeff, t(k), u(k,:,:))
       deallocate(lincoeff)        
    end if

  end subroutine step

  !================================================================!
  ! Retrieve the coefficients for linearizing the jacobian
  !================================================================!
  
  impure subroutine get_linearization_coeff(this, cindex, h, lincoeff)

    class(ABM)    , intent(in)    :: this
    type(integer) , intent(in)    :: cindex
    type(scalar)  , intent(in)    :: h 
    type(scalar)  , intent(inout) :: lincoeff(:)
    type(integer)                 :: p
    
    associate(&
         & deriv_order => this % system % get_differential_order(), &
         & a => this % A(cindex,1))
      
      forall(p = 0:deriv_order)
         lincoeff(p+1) = (a*h)**(deriv_order-p)
      end forall
      
    end associate

  end subroutine get_linearization_coeff

end module abm_integrator_class
