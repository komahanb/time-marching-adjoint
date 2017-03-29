!=====================================================================!
! Main Program for testing the integrators on different test problems
!=====================================================================!

#include "scalar.fpp"

program test_time_integration

  use vanderpol_system          , only : fvanderpol => vanderpol_first_order 
  use dynamic_physics_interface , only : dynamics

  implicit none

  class(dynamics), allocatable :: sys

  test_vanderpol: block
    allocate(sys, source = fvanderpol(0.0d0))
    call test_integrators(sys)
    deallocate(sys)
  end block test_vanderpol

contains

  subroutine test_integrators( test_system)

    use abm_integrator_class     , only : ABM

    class(dynamics), intent(inout) :: test_system    
    type(ABM)                      :: abmobj

    ! Create the integrator
    abmobj = ABM(system = test_system, tinit=0.0d0, tfinal = 10.0d0, &
         & h=1.0d-3, implicit=.true., max_abm_order=1)

    call abmobj % set_print_level(2)
    call abmobj % to_string()
    call abmobj % integrate()
    call abmobj % write_solution("abm.dat")
    call abmobj % to_string()

  end subroutine test_integrators

end program test_time_integration

