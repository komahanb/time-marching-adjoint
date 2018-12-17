!=====================================================================!
! Main Program for testing the integrators on different test problems
!=====================================================================!

#include "scalar.fpp"

program test

  use dynamic_physics_interface             , only : dynamics
  use backward_differences_integrator_class , only : bdf
  use decay_ode_class                       , only : decay

  implicit none

  class(dynamics), allocatable :: system
  type(bdf)                    :: bdfobj

  allocate(system, source = decay(gamma = 0.0d0))
  bdfobj = BDF(system = system, tinit=0.0d0, tfinal = 1.0d0, &
       & h=1.0d-3, implicit=.true., accuracy_order=2)
  call bdfobj % to_string()
  call bdfobj % solve()
  call bdfobj % write_solution("decay-bdf.dat")
  deallocate(system)

contains
  
end program test
