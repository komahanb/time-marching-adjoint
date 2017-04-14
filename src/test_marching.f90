!=====================================================================!
! Main Program for testing the integrators on different test problems
!=====================================================================!

#include "scalar.fpp"

program test_time_integration

  use vanderpol_system          , only : fvanderpol => vanderpol_first_order
  use spring_dynamics_class     , only : smd
  use freefall_dynamics_class   , only : freefall
  use test_ode_class            , only : ODE
  use dynamic_physics_interface , only : dynamics

  implicit none

  class(dynamics), allocatable :: sys

  test_vanderpol: block
    !allocate(sys, source = smd(2.0d0, 0.0d0, 2.0d0))
    !allocate(sys, source = fvanderpol(1.0d0))
    !allocate(sys, source = freefall(1.0d0, -10.0d0))
    allocate(sys, source = ODE(A=[2.0d0, 2.0d0, 2.0d0], order=2, nvars=3))
    call sys % set_approximate_jacobian(.false.)    
    call test_integrators(sys)
    deallocate(sys)
  end block test_vanderpol

contains

  subroutine test_integrators( test_system)

    use abm_integrator_class     , only : ABM
    use newmark_integrator_class , only : newmark
    use runge_kutta_integrator_class , only : dirk
    use backward_differences_integrator_class , only : bdf

    class(dynamics), intent(inout) :: test_system    
    type(ABM)     :: abmobj
    type(newmark) :: nbg
    type(dirk)    :: dirkobj
    type(bdf)     :: bdfobj

    abmobj = ABM(system = test_system, tinit=0.0d0, tfinal = 10.0d0, &
         & h=1.0d-3, implicit=.true., max_order=3)
    call abmobj % to_string()
    call abmobj % integrate()
    call abmobj % write_solution("abm.dat")
    call abmobj % to_string()

    nbg = newmark(system = test_system, tinit=0.0d0, tfinal = 10.0d0, &
         & h=1.0d-3, implicit=.true., max_order=3)
    call nbg % to_string()
    call nbg % integrate()
    call nbg % write_solution("nbg.dat")
    call nbg % to_string()

    bdfobj = BDF(system = test_system, tinit=0.0d0, tfinal = 10.0d0, &
         & h=1.0d-3, implicit=.true., max_order=6)
    call bdfobj % to_string()
    call bdfobj % integrate()
    call bdfobj % write_solution("bdf.dat")
    call bdfobj % to_string()

    dirkobj = DIRK(system = test_system, tinit=0.0d0, tfinal = 10.0d0, &
         & h=1.0d-3, implicit=.true., max_order=2)
    call dirkobj % to_string()
    call dirkobj % integrate()
    call dirkobj % write_solution("dirk.dat")
    call dirkobj % to_string()
    
  end subroutine test_integrators

end program test_time_integration

