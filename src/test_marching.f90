!=====================================================================!
! Main Program for testing the integrators on different test problems
!=====================================================================!

#include "scalar.fpp"

program test_time_integration
  
  ! Import all test physics
  !use spring_mass_damper_class      , only : smd1, smd2
  use vanderpol_system , only : fvanderpol => vanderpol_first_order 
  use dynamic_physics_interface, only : dynamics

  implicit none
  
  class(dynamics), allocatable :: sys

!!$  ! Declare Physics for testing
!!$  type(smd1)      , target :: smd1obj                 ! Spring-mass-damper test ODE (1 var)
!!$  type(smd2)      , target :: smd2obj                 ! Spring-mass-damper test ODE (2 var)
!!$  type(vanderpol) , target :: vpl                     ! Vanderpol equation (2 var)
!!$  type(aero_elastic_oscillator), target :: aeosc      ! Aeroelastic oscillator (2 vars)

  ! Test the integrators in vanderpol oscillator system  
  test_vanderpol: block   
    allocate(sys, source = fvanderpol(1.0d0))
    call test_integrators(sys)
  end block test_vanderpol
  
!!$  ! Test the integrators on 2 dof spring mass damper system
!!$  test_smd2: block
!!$    call smd2obj % initialize("SMD2", num_state_vars = 2)
!!$    call test_integrators(smd2obj, "smd2", .true.)
!!$    call smd2obj % finalize()
!!$  end block test_smd2
!!$
!!$  ! Test the integrators on aeroelastic oscillator problem with
!!$  ! pitching and plunging degree of freedom
!!$  test_aeosc: block
!!$    call aeosc % initialize("AeroElasticOscillator", num_state_vars = 2)
!!$    call test_integrators(aeosc, "aeosc", .true.)
!!$    call aeosc % finalize()
!!$  end block test_aeosc

contains

  subroutine test_integrators( test_system)

    ! Import all time integrators
!!$    use runge_kutta_integrator        , only : DIRK
!!$    use bdf_integrator                , only : BDF
    use abm_integrator_class     , only : ABM
!!$    use nbg_integrator                , only : NBG
    
    class(dynamics), intent(inout) :: test_system    

    ! Declare integrators
    ! type(DIRK) :: dirkobj   ! DIRK Integrator object
    ! type(BDF)  :: bdfobj    ! BDF Integrator object
    type(ABM)  :: abmobj    ! ABM Integrator object
    !type(NBG)  :: nbgobj    ! NBM Integrator object
!!$    
!!$    !=================================================================!
!!$    !                        TEST BDF                                 !
!!$    !=================================================================!
!!$
!!$    bdfobj = BDF(system = test_system, tfinal = 20.0d0, h=1.0d-2, &
!!$         & max_bdf_order = 3, second_order=second_order)
!!$    call bdfobj % set_print_level(2)
!!$    call bdfobj % integrate()
!!$    call bdfobj % write_solution("_bdf.dat")
!!$
!!$    !=================================================================!
!!$    !                     TEST NBG                                    !
!!$    !=================================================================!
!!$
!!$    nbgobj = NBG(system = test_system, tfinal = 20.0d0, h=1.0d-2, &
!!$         & second_order=second_order)
!!$    call nbgobj % set_print_level(2)
!!$    call nbgobj % integrate()
!!$    call nbgobj % write_solution("_nbg.dat")

    !=================================================================!
    !                     TEST ABM                                    !
    !==================================================================!

    abmobj = ABM(system = test_system, tinit=0.0d0, tfinal = 1.0d0, h=1.0d-1, implicit=.true., max_abm_order=2)
    call abmobj % set_print_level(2)
    call abmobj % to_string()
    call abmobj % integrate()
    call abmobj % write_solution("_abm.dat")

!!$
!!$    !=================================================================!
!!$    !                        TEST DIRK                                !
!!$    !=================================================================!
!!$
!!$    dirkobj = DIRK(system = test_system, tfinal = 20.0d0, h=1.0d-2, &
!!$         & num_stages=2, second_order=second_order) 
!!$    call dirkobj % set_print_level(2)
!!$    call dirkobj % integrate()
!!$    call dirkobj % write_solution(name//"_dirk.dat")

  end subroutine test_integrators
  
end program test_time_integration

