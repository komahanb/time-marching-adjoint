!=====================================================================!
! Main Program for testing the integrators on different test problems
!=====================================================================!

#include "scalar.fpp"

program test_time_integration

  use vanderpol_system          , only : fvanderpol => vanderpol_first_order
  use spring_dynamics_class     , only : smd
  use freefall_dynamics_class   , only : freefall
  use test_ode_class            , only : ODE
  use first_order_ode_class     , only : FODE
  use dynamic_physics_interface , only : dynamics

  implicit none

  class(dynamics), allocatable :: sys

  test_vanderpol: block
    !allocate(sys, source = smd(0.5d0, 6.0d0, 18.0d0))
    allocate(sys, source = FODE(9.8d0, -0.196d0))
    !allocate(sys, source = fvanderpol(1.0d0))
    !allocate(sys, source = freefall(1.0d0, -10.0d0))
    !allocate(sys, source = ODE(A=[2.0d0, 2.0d0, 2.0d0], order=4, nvars=3))
    call test_integrators(sys, 0.0d0, 10.0d0, 1250, "-1250")
    call test_integrators(sys, 0.0d0, 10.0d0, 2500, "-2500")
    call test_integrators(sys, 0.0d0, 10.0d0, 5000, "-5000")
    call test_integrators(sys, 0.0d0, 10.0d0, 10000, "-10000")
    deallocate(sys)
  end block test_vanderpol

contains

  subroutine test_integrators(test_system, tinit, tfinal, nsteps, suffix)

    use abm_integrator_class , only : ABM
    use newmark_integrator_class , only : newmark
    use runge_kutta_integrator_class , only : dirk
    use backward_differences_integrator_class , only : bdf

    class(dynamics), intent(inout) :: test_system
    type(scalar)  :: tinit, tfinal
    type(integer) :: nsteps
    type(ABM)     :: abmobj
    type(newmark) :: nbg
    type(dirk)    :: dirkobj
    type(bdf)     :: bdfobj       
    type(scalar)  :: h
    character(len=*) :: suffix

    h = (tfinal-tinit)/dble(nsteps)

    ! Dirk
    dirkobj = DIRK(system = test_system, tinit=tinit, tfinal = tfinal, &
         & h=h, implicit=.true., accuracy_order=2)
    call dirkobj % to_string()
    call dirkobj % solve()
    call dirkobj % write_solution("smd-dirk2"//suffix//".dat")

    dirkobj = DIRK(system = test_system, tinit=tinit, tfinal = tfinal, &
         & h=h, implicit=.true., accuracy_order=3)
    call dirkobj % to_string()
    call dirkobj % solve()
    call dirkobj % write_solution("smd-dirk3"//suffix//".dat")

    dirkobj = DIRK(system = test_system, tinit=tinit, tfinal = tfinal, &
         & h=h, implicit=.true., accuracy_order=4)
    call dirkobj % to_string()
    call dirkobj % solve()
    call dirkobj % write_solution("smd-dirk4"//suffix//".dat")

    ! bdf
    bdfobj = BDF(system = test_system, tinit=tinit, tfinal = tfinal, &
         & h=h, implicit=.true., accuracy_order=1)
    call bdfobj % to_string()
    call bdfobj % solve()
    call bdfobj % write_solution("smd-bdf1"//suffix//".dat")

    bdfobj = BDF(system = test_system, tinit=tinit, tfinal = tfinal, &
         & h=h, implicit=.true., accuracy_order=2)
    call bdfobj % to_string()
    call bdfobj % solve()
    call bdfobj % write_solution("smd-bdf2"//suffix//".dat")

    bdfobj = BDF(system = test_system, tinit=tinit, tfinal = tfinal, &
         & h=h, implicit=.true., accuracy_order=3)
    call bdfobj % to_string()
    call bdfobj % solve()
    call bdfobj % write_solution("smd-bdf3"//suffix//".dat")

    bdfobj = BDF(system = test_system, tinit=tinit, tfinal = tfinal, &
         & h=h, implicit=.true., accuracy_order=4)
    call bdfobj % to_string()
    call bdfobj % solve()
    call bdfobj % write_solution("smd-bdf4"//suffix//".dat")

    bdfobj = BDF(system = test_system, tinit=tinit, tfinal = tfinal, &
         & h=h, implicit=.true., accuracy_order=5)
    call bdfobj % to_string()
    call bdfobj % solve()
    call bdfobj % write_solution("smd-bdf5"//suffix//".dat")

    bdfobj = BDF(system = test_system, tinit=tinit, tfinal = tfinal, &
         & h=h, implicit=.true., accuracy_order=6)
    call bdfobj % to_string()
    call bdfobj % solve()
    call bdfobj % write_solution("smd-bdf6"//suffix//".dat")

    ! abm
    abmobj = ABM(system = test_system, tinit=tinit, tfinal = tfinal, &
         & h=h, implicit=.true., accuracy_order=1)
    call abmobj % to_string()
    call abmobj % solve()
    call abmobj % write_solution("smd-abm1"//suffix//".dat")

    abmobj = ABM(system = test_system, tinit=tinit, tfinal = tfinal, &
         & h=h, implicit=.true., accuracy_order=2)
    call abmobj % to_string()
    call abmobj % solve()
    call abmobj % write_solution("smd-abm2"//suffix//".dat")

    abmobj = ABM(system = test_system, tinit=tinit, tfinal = tfinal, &
         & h=h, implicit=.true., accuracy_order=3)
    call abmobj % to_string()
    call abmobj % solve()
    call abmobj % write_solution("smd-abm3"//suffix//".dat")

    abmobj = ABM(system = test_system, tinit=tinit, tfinal = tfinal, &
         & h=h, implicit=.true., accuracy_order=4)
    call abmobj % to_string()
    call abmobj % solve()
    call abmobj % write_solution("smd-abm4"//suffix//".dat")

    abmobj = ABM(system = test_system, tinit=tinit, tfinal = tfinal, &
         & h=h, implicit=.true., accuracy_order=5)
    call abmobj % to_string()
    call abmobj % solve()
    call abmobj % write_solution("smd-abm5"//suffix//".dat")

    abmobj = ABM(system = test_system, tinit=tinit, tfinal = tfinal, &
         & h=h, implicit=.true., accuracy_order=6)
    call abmobj % to_string()
    call abmobj % solve()
    call abmobj % write_solution("smd-abm6"//suffix//".dat")

  end subroutine test_integrators

end program test_time_integration

