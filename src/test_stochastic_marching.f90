!=====================================================================!
! Main Program for testing the integrators on different test problems
!=====================================================================!

#include "scalar.fpp"

program test

  use iso_fortran_env                       , only : dp => REAL64
  use dynamic_physics_interface             , only : dynamics
  use backward_differences_integrator_class , only : bdf
  use decay_ode_class                       , only : decay
  use class_random_normal                   , only : normal

  implicit none

  integer, parameter :: nmcs = 1000
  integer            :: imcs
  character*60       :: cmcs

  class(dynamics), allocatable :: system
  type(bdf)    :: bdfobj
  type(normal) :: nrand
  real(dp)     :: gamma

  real(dp), allocatable :: umcs(:,:)
  real(dp), allocatable :: umean(:)
  real(dp), allocatable :: uvar(:)

  ! Initialize random normal number generator
  nrand = normal(mu=0.0d0, sigma=1.0d-1)
  
  allocate(umcs(1001,nmcs)) 
  umcs = 0

  do imcs = 1, nmcs

     write(cmcs,*) imcs
     call nrand % get(gamma)
     allocate(system, source = decay(gamma))
     bdfobj = BDF( system = system, tinit=0.0d0, tfinal = 1.0d0, &
          & h=1.0d-3, implicit=.true., accuracy_order=2 )
     print *, imcs, gamma
     call bdfobj % solve()
     umcs(:,imcs) = bdfobj % u(:,1,1)
     !call bdfobj % write_solution("mcs-"//trim(adjustl(cmcs))//".dat")
     deallocate(system)

  end do

  allocate(umean(1001))
  umean = 0
  allocate(uvar(1001))
  uvar = 0
 
  allocate(system, source = decay(gamma))
  umean = sum(umcs, dim=2)/dble(nmcs)
  bdfobj % u (:,1,1) = umean
  call bdfobj % write_solution("mcs-mean.dat")
  
end program test
