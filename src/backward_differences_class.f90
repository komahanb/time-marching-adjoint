#include "scalar.fpp"

!=====================================================================!
! Backward Difference Formulas integration module for differential 
! equations.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================! 

module backward_differences_integrator_class

  use integrator_interface      , only : integrator
  use dynamic_physics_interface , only : dynamics
  use utils, only : real_part

  implicit none

  private
  public :: BDF
  
  !===================================================================! 
  ! BDF Integrator type
  !===================================================================! 

  type, extends(integrator) :: BDF
     
     private

     ! BDF variables
     type(integer)             :: max_bdf_order = 6
     type(scalar), allocatable :: A(:,:)

   contains
           
     procedure :: evaluate_states
     procedure :: get_linear_coeff
     procedure :: get_order

     ! Destructor
     final :: destroy

  end type BDF

  interface BDF
     module procedure create
  end interface BDF

contains

  !===================================================================!
  ! Initialize the BDF datatype and allocate required variables
  !===================================================================!
  
  type(bdf) function create(system, tinit, tfinal, h, implicit, &
       & max_order) result(this)

    class(dynamics)   , intent(in)   , target :: system
    type(scalar)      , intent(in)            :: tinit, tfinal
    type(scalar)      , intent(in)            :: h
    type(integer)     , intent(in)            :: max_order
    type(logical)     , intent(in)            :: implicit   

    print *, "======================================"
    print *, ">>   Backward Difference Formulas  << "
    print *, "======================================"

    call this % construct(system, tinit, tfinal, h, implicit, 0)

    !-----------------------------------------------------------------!
    ! Set the order of integration
    !-----------------------------------------------------------------!

    if (max_order .le. this % max_bdf_order) this % max_bdf_order = max_order
    print '("  >> Max BDF Order          : ",i4)', this % max_bdf_order

    allocate( this % A (this % max_bdf_order, this % max_bdf_order+1) )
    this % A = 0.0d0 

    ! Set the coefficients
    ! http://www.scholarpedia.org/article/Backward_differentiation_formulas
    if ( this % max_bdf_order .ge. 1 ) this % A(1,1:2) = [1.0d0, -1.0d0]
    if ( this % max_bdf_order .ge. 2 ) this % A(2,1:3) = [3.0d0, -4.0d0, 1.0d0]/2.0d0
    if ( this % max_bdf_order .ge. 3 ) this % A(3,1:4) = [11.0d0, -18.0d0, 9.0d0, -2.0d0]/6.0d0
    if ( this % max_bdf_order .ge. 4 ) this % A(4,1:5) = [25.0d0, -48.0d0, 36.0d0, -16.0d0, 3.0d0]/12.0d0
    if ( this % max_bdf_order .ge. 5 ) this % A(5,1:6) = [137.0d0, -300.0d0, 300.0d0, -200.0d0, 75.0d0, -12.0d0]/60.0d0
    if ( this % max_bdf_order .ge. 6 ) this % A(6,1:7) = [147.0d0, -360.0d0, 450.0d0, -400.0d0, 225.0d0, -72.0d0, 10.0d0]/60.0d0

    ! Sanity check on BDF coeffs
    sanity_check: block
      type(integer) :: j
      do j = 1, this % max_bdf_order
         if (abs(sum(real_part(this % A(j,1:j+1)))) .gt. 1.0d-15 ) then
            print *, "Error in BDF Coeff for order ", abs(sum(real_part(this % A(j,1:j+1)))), j
            stop
         end if
      end do
    end block sanity_check
    
  end function create

  !=================================================================!
  ! Destructor for the BDF integrator
  !=================================================================!
  
  impure subroutine destroy(this)

    type(BDF), intent(inout) :: this

    ! Parent class call
    call this % destruct()

    ! Deallocate BDF coefficient
    if(allocated(this % A)) deallocate(this % A)

  end subroutine destroy

  !===================================================================!
  ! Returns the order of approximation for the given time step k and
  ! degree d
  !===================================================================!

  impure type(integer) function get_order(this, step) result(order)

    class(BDF)   , intent(in) :: this
    type(integer), intent(in) :: step

    order = step - 1

    if (order .gt. this % max_bdf_order) order = this % max_bdf_order

  end function get_order

  !================================================================!
  ! Interface to approximate states using the time marching coeffs
  !================================================================!

  impure subroutine evaluate_states(this, t, u)

    use nonlinear_algebra, only : nonlinear_solve

    class(BDF)   , intent(in)    :: this
    type(scalar) , intent(in)    :: t(:)      ! array of time values
    type(scalar) , intent(inout) :: u(:,:,:)  ! previous values of state variables

    type(scalar)  , allocatable :: lincoeff(:)  ! order of equation + 1
    type(integer) :: k , i, n
    type(scalar)  :: scale
    
    ! Pull out the number of time steps of states provided and add one
    ! to point to the current time step
    k = size(u(:,1,1))

    associate( &
         & p => this % get_order(k), &
         & A => this % A(this % get_order(k),:), &
         & h => t(k) - t(k-1), &
         & torder => this % system % get_time_deriv_order())

      ! Assume a value for lowest order state
      u(k,1,:) = 0

      ! Find the higher order states based on BDF formula
      do n = 1, torder
         do i = 0, p
            scale = A(i+1)/h !(t(k-i)-t(k-i-1))
            u(k,n+1,:) = u(k,n+1,:) + scale*u(k-i,n,:)
         end do
      end do

      ! Perform a nonlinear solution if this is a implicit method
      if ( this % is_implicit() ) then
         allocate(lincoeff(torder + 1))
         call this % get_linear_coeff(lincoeff, p, h)
         call nonlinear_solve(this % system, lincoeff, t(k), u(k,:,:))
         deallocate(lincoeff)
      end if

    end associate

end subroutine evaluate_states

  !================================================================!
  ! Retrieve the coefficients for linearizing the jacobian
  !================================================================!
  
  impure subroutine get_linear_coeff(this, lincoeff, int_order, h)

    class(BDF)    , intent(in)  :: this
    type(integer) , intent(in)  :: int_order     ! order of approximation of the integration
    type(scalar)  , intent(in)  :: h ! step size
    type(scalar)  , intent(inout) :: lincoeff(:)   ! order of equation + 1   
    type(integer) :: p
    
    associate(&
         & deriv_order => this % system % get_time_deriv_order(), &
         & a => this % A(int_order,1))
      
      forall(p = 0:deriv_order)
         lincoeff(p+1) = (a/h)**p
      end forall
      
    end associate

  end subroutine get_linear_coeff

end module backward_differences_integrator_class
