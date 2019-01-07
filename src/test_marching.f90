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
    allocate(sys, source = smd(1.0d0, 0.5d0, 5.0d0))
    !allocate(sys, source = fvanderpol(1.0d0))
    !allocate(sys, source = freefall(1.0d0, -10.0d0))
    !allocate(sys, source = ODE(A=[2.0d0, 2.0d0, 2.0d0], order=4, nvars=3))
    call test_integrators(sys)
    deallocate(sys)
  end block test_vanderpol

contains

  subroutine test_integrators(test_system)

    use abm_integrator_class , only : ABM
    use newmark_integrator_class , only : newmark
    use runge_kutta_integrator_class , only : dirk
    use backward_differences_integrator_class , only : bdf

    class(dynamics), intent(inout) :: test_system    
    type(ABM)     :: abmobj
    type(newmark) :: nbg
    type(dirk)    :: dirkobj
    type(bdf)     :: bdfobj
        
!!$    abmobj = ABM(system = test_system, tinit=0.0d0, tfinal = 10.0d0, &
!!$         & h=1.0d-3, implicit=.true., accuracy_order=6)
!!$    call abmobj % to_string()
!!$    call abmobj % solve()
!!$    call abmobj % write_solution("abm.dat")

    dirkobj = DIRK(system = test_system, tinit=0.0d0, tfinal = 10.0d0, &
         & h=1.0d-3, implicit=.true., accuracy_order=2)
    call dirkobj % to_string()
    call dirkobj % solve()
    call dirkobj % write_solution("smd-dirk2.dat")

    dirkobj = DIRK(system = test_system, tinit=0.0d0, tfinal = 10.0d0, &
         & h=1.0d-3, implicit=.true., accuracy_order=3)
    call dirkobj % to_string()
    call dirkobj % solve()
    call dirkobj % write_solution("smd-dirk3.dat")

    dirkobj = DIRK(system = test_system, tinit=0.0d0, tfinal = 10.0d0, &
         & h=1.0d-3, implicit=.true., accuracy_order=4)
    call dirkobj % to_string()
    call dirkobj % solve()
    call dirkobj % write_solution("smd-dirk4.dat")

!!$    bdfobj = BDF(system = test_system, tinit=0.0d0, tfinal = 10.0d0, &
!!$         & h=1.0d-3, implicit=.true., accuracy_order=6)
!!$    call bdfobj % to_string()
!!$    call bdfobj % solve()
!!$    call bdfobj % write_solution("bdf.dat")
    
!!$    if ( test_system % get_differential_order() .eq. 2 ) then
!!$       nbg = newmark(system = test_system, tinit=0.0d0, tfinal = 10.0d0, &
!!$            & h=1.0d-3, implicit=.true., accuracy_order=2)
!!$       call nbg % to_string()
!!$       call nbg % solve()
!!$       call nbg % write_solution("nbg.dat")
!!$    end if    

  end subroutine test_integrators

end program test_time_integration

