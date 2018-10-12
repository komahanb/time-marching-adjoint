!=====================================================================!
! Solve one dimensional transport equations using integrators
!=====================================================================!

#include "scalar.fpp"

program test_time_integration

  use constants                 , only : WP
  use dynamic_physics_interface , only : dynamics
  use unsteady_transport_class  , only : unsteady_transport

  implicit none

  class(dynamics), allocatable :: system
  type(scalar)   , parameter   :: bounds(2) = [5.0_wp, 45.0_wp]

  test_transport: block
    allocate(system, source = unsteady_transport( &
         & diffusion_coeff = 0.01_WP, &
         & convective_velocity = 1.0_WP, &
         & bounds = bounds, npts=100 ))
    call test_integrators(system)
    deallocate(system)
  end block test_transport

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
        
    abmobj = ABM(system = test_system, tinit=0.0d0, tfinal = 10.0d0, &
         & h=1.0d-3, implicit=.true., accuracy_order=6)
    call abmobj % to_string()
    call abmobj % solve()
    call abmobj % write_solution("transport-abm.dat")

    dirkobj = DIRK(system = test_system, tinit=0.0d0, tfinal = 10.0d0, &
         & h=1.0d-3, implicit=.true., accuracy_order=4)
    call dirkobj % to_string()
    call dirkobj % solve()
    call dirkobj % write_solution("transport-dirk.dat")
    
    bdfobj = BDF(system = test_system, tinit=0.0d0, tfinal = 10.0d0, &
         & h=1.0d-3, implicit=.true., accuracy_order=6)
    call bdfobj % to_string()
    call bdfobj % solve()
    call bdfobj % write_solution("transport-bdf.dat")
    
    if ( test_system % get_differential_order() .eq. 2 ) then
       nbg = newmark(system = test_system, tinit=0.0d0, tfinal = 10.0d0, &
            & h=1.0d-3, implicit=.true., accuracy_order=2)
       call nbg % to_string()
       call nbg % solve()
       call nbg % write_solution("transport-nbg.dat")
    end if    

  end subroutine test_integrators

end program test_time_integration

