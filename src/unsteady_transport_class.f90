#include "scalar.fpp"

!=====================================================================!
! One dimensional unsteady tranport physics. 
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module unsteady_transport_class

  use constants                 , only : WP
  use dynamic_physics_interface , only : dynamics

  implicit none

  private

  public :: unsteady_transport

  !-------------------------------------------------------------------!
  ! Type that implements first order transport equations
  !-------------------------------------------------------------------!
  
  type, extends(dynamics) :: unsteady_transport

     type(scalar) :: dx
     type(scalar) :: conv_speed
     type(scalar) :: diff_coeff 
     
   contains
     
     ! Implement deferred procedures from superclasses
     procedure :: add_residual
     procedure :: add_jacobian
     procedure :: get_initial_condition
     
     ! Destructor
     final :: destruct
     
  end type unsteady_transport
  
  ! Interface to construct the physical system
  interface unsteady_transport
     procedure construct_unsteady_transport
  end interface unsteady_transport

contains
 
  !===================================================================!
  ! Constructor for unsteady transport physics
  !===================================================================!
  
  type(unsteady_transport) function construct_unsteady_transport( &
       & diffusion_coeff, convective_velocity, &
       & bounds, npts, sparse) result (this)

    type(scalar), intent(in) :: diffusion_coeff, convective_velocity
    type(scalar), intent(in) :: bounds(2)
    integer     , intent(in) :: npts
    logical     , intent(in) :: sparse

    type(scalar) :: dx
    integer      :: i

    call this % set_description('Unsteady Transport')
    
    ! Set time order of physical system
    call this % set_differential_order(1)

    ! System parameters
    this % conv_speed = convective_velocity
    this % diff_coeff = diffusion_coeff

    ! Sparse or dense assembly
    this % sparse = sparse

    allocate(this % X(3,npts+2)); this % X = 0.0_wp;
    this % dx = (bounds(2) - bounds(1))/dble(npts+1)
    do i = 1, npts + 2
       this % X(1,i) = bounds(1) + dble(i-1)*this % dx
    end do
    
    ! Set the number of state variables based on spatial
    ! discretization of the governing equations
    call this % set_num_state_vars(npts) ! ignore boundary nodes
    
  end function construct_unsteady_transport
  
  !===================================================================!
  ! Destructor for unsteady transport physics
  !===================================================================!
  
  pure subroutine destruct(this)
    
    type(unsteady_transport), intent(inout) :: this
    
  end subroutine destruct
  
  !===================================================================!
  ! Residual assembly at each time step
  !===================================================================!
  
  pure subroutine add_residual(this, residual, U, X)

    class(unsteady_transport), intent(inout) :: this
    type(scalar)             , intent(inout) :: residual(:)
    type(scalar)             , intent(in)    :: U(:,:)
    type(scalar)             , intent(in)    :: X(:,:)    

    ! Locals
    type(scalar) :: a, b, c
    integer :: i

    associate(phi=>U(1,:), phidot=> U(2,:), &
         & npts=>this % get_num_state_vars(), &
         & vel=>this % conv_speed, gamma => this % diff_coeff)

    a = -(vel/(2.0_wp*this % dx) + gamma/(this % dx*this % dx))
    b = 2.0_wp*gamma/(this % dx*this % dx)
    c = (vel/(2.0_wp*this % dx) - gamma/(this % dx*this % dx))

    ! Residual with BC applied
    residual(1) = residual(1) +  b*phi(1) + c*phi(2) + phidot(1) -a*0.0_wp
    forall(i = 2:npts-1)
       residual(i) = residual(i) + a*phi(i-1) + b*phi(i) + c*phi(i+1) + phidot(i)
    end forall
    residual(npts) = residual(npts) + a*phi(npts-1) + b*phi(npts) + phidot(npts) -c*0.0_wp

  end associate

end subroutine add_residual

  !===================================================================!
  ! Jacobian assembly at each time step. 
  !===================================================================!
  
  pure subroutine add_jacobian(this, jacobian, coeff, U, X)

    class(unsteady_transport) , intent(inout) :: this
    type(scalar)              , intent(inout) :: jacobian(:,:)
    type(scalar)              , intent(in)    :: coeff(:)
    type(scalar)              , intent(in)    :: U(:,:)
    type(scalar)              , intent(in)    :: X(:,:)    

    ! Locals
    type(scalar) :: aa, bb, cc
    integer      :: i, j
    
    associate(phi=>U(1,:), phidot=> U(2,:), &
         & npts=>this % get_num_state_vars(), &
         & vel=>this % conv_speed, gamma => this % diff_coeff, &
         & alpha=>coeff(1), beta=>coeff(2))      

      aa = -(vel/(2.0_wp*this % dx) + gamma/(this % dx*this % dx))
      bb = 2.0_wp*gamma/(this % dx*this % dx)
      cc = (vel/(2.0_wp*this % dx) - gamma/(this % dx*this % dx))

      if (this % sparse .eqv. .true.) then           

         jacobian(1,:) = [0.0d0, beta + alpha*bb, alpha*cc]
         do concurrent(i = 2 : npts-1)
            jacobian(i,:) = [alpha*aa, beta + alpha*bb, alpha*cc]
         end do
         jacobian(npts,:) = [alpha*aa, beta + alpha*bb, 0.0d0]

      else

         do i = 1, npts
            do j = 1, npts
               if (i .eq. j-1) then
                  ! upper diagonal
                  jacobian(i,j) = jacobian(i,j) + alpha*cc 
               else if (i .eq. j) then
                  ! diagonal
                  jacobian(i,i) = jacobian(i,j) + beta + alpha*bb
               else if (i .eq. j+1) then
                  ! lower diagonal
                  jacobian(i,j) = jacobian(i,j) + alpha*aa
               else
               end if
            end do
         end do

      end if

    end associate

  end subroutine add_jacobian
  
  !===================================================================!
  ! Sets the initial condition for use in the integator. 
  !===================================================================!  

  pure subroutine get_initial_condition(this, U, X)
    
    class(unsteady_transport), intent(in)    :: this
    type(scalar)             , intent(inout) :: U(:,:)
    type(scalar)             , intent(in)    :: X(:,:)
    type(scalar) :: pi
    integer :: i

    pi = 4.0_wp*atan(1.0_wp)
    
    associate(phi=>U(1,:), x=>this % X(1,:))
      
      do i=1,this % get_num_state_vars()
         phi(i) = (0.4_wp*pi)**(-0.5_wp)*exp(-2.5_wp*(x(i)-10.0_wp)*(x(i)-10.0_wp))
      end do
      
    end associate
    
  end subroutine get_initial_condition

end module unsteady_transport_class
