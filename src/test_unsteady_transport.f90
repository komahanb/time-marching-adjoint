!=====================================================================!
! Solve one dimensional transport equations using integrators
!=====================================================================!

#include "scalar.fpp"

program test_time_integration

  use constants                 , only : WP
  use dynamic_physics_interface , only : dynamics
  use unsteady_transport_class  , only : unsteady_transport
  
  implicit none

  test_transport: block

  class(dynamics), allocatable :: system
  type(scalar)   , parameter   :: bounds(2) = [5.0_wp, 45.0_wp]

    ! case 1
    allocate(system, source = unsteady_transport( &
         & diffusion_coeff = 0.01_WP, &
         & convective_velocity = 1.0_WP, &
         & bounds = bounds, npts=500, &
         & sparse = .true.))
    call test_integrators(system, 'case1')
    deallocate(system)
    
    !case 2
    allocate(system, source = unsteady_transport( &
         & diffusion_coeff = 0.0_WP, &
         & convective_velocity = 1.0_WP, &
         & bounds = bounds, npts=500, &
         & sparse = .true.))
    call test_integrators(system, 'case2')
    deallocate(system)

  end block test_transport

  spatial_convergence : block

    real(wp) :: rmse(6), h
    integer  :: npts = 100
    integer  :: n

    ! grid spacing
    open(12, file='spatial-error.dat')
    write(12, *) "npts ", "h ", "dirk2 ", "dirk3 ", "dirk4 ", "exp-euler ", "imp-euler ", "cni "

    do n = 1, 4
       h = 40.0d0/dble(npts+1)
       call evaluate_spatial_error(npts, rmse)
       write(12, *) npts, h, rmse(1), rmse(2), rmse(3), rmse(4), rmse(5), rmse(6)
       npts = npts*10
    end do
    close(12)

  end block spatial_convergence

!!$  temporal_convergence : block
!!$
!!$    real(wp) :: rmse(3), h
!!$    integer  :: nsteps = 1000
!!$    integer  :: n
!!$
!!$    ! grid spacing
!!$    open(12, file='imp-euler-temporal-error.dat')
!!$    write(12, *) "nsteps ", "h ", "rmse1 ", "rmse2 ", "rmse3 "
!!$
!!$    do n = 1, 5
!!$       h = 5.0d0/dble(nsteps+1)
!!$       call evaluate_temporal_error(nsteps, rmse)
!!$       write(12, *) nsteps, h, rmse(1), rmse(2), rmse(3)
!!$       nsteps = nsteps*2
!!$    end do
!!$
!!$    close(12)
!!$
!!$  end block temporal_convergence

contains

!!$  subroutine evaluate_temporal_error(nsteps, rmse)
!!$
!!$    use dynamic_physics_interface             , only : dynamics
!!$    use backward_differences_integrator_class , only : bdf
!!$    use unsteady_transport_class              , only : unsteady_transport   
!!$
!!$    integer         , intent(in)  :: nsteps
!!$    real(wp)        , intent(inout) :: rmse(3)
!!$    type(scalar)    , parameter   :: bounds(2) = [5.0_wp, 45.0_wp]
!!$    class(dynamics) , allocatable :: system
!!$    type(bdf)                     :: bdfobj
!!$    integer                       :: j, k
!!$    real(wp)                      :: error, t, x
!!$    real(wp) :: h
!!$
!!$    h = (15.0_wp - 10.0_wp)/dble(nsteps+1)
!!$
!!$    ! Case 1
!!$    allocate(system, source = unsteady_transport( &
!!$         & diffusion_coeff = 0.01_WP, &
!!$         & convective_velocity = 1.0_WP, &
!!$         & bounds = bounds, npts = 1000, &
!!$         & sparse = .true.))
!!$
!!$    ! create integrator with supplied mesh resolution
!!$    bdfobj = BDF(system = system, tinit=10.0d0, tfinal = 15.0d0, &
!!$         & h=h, implicit = .false., accuracy_order=1)
!!$    call bdfobj % to_string()
!!$    call bdfobj % solve()
!!$
!!$    ! find rmse at the last time step
!!$    rmse(1) = 0.0d0   
!!$    do j = 1, system % get_num_state_vars()
!!$       do k = 1, bdfobj % num_time_steps
!!$          t = bdfobj % time(k)
!!$          x = system % X(1,j+1)
!!$          error = bdfobj % U (k, 1, j) - exact(t, x)
!!$          rmse(1) = rmse(1) + error**2.0_wp
!!$       end do
!!$    end do
!!$    rmse(1) = sqrt(rmse(1)/(system % get_num_state_vars()*bdfobj % num_time_steps))
!!$
!!$    ! create integrator with supplied mesh resolution
!!$    bdfobj = BDF(system = system, tinit=10.0d0, tfinal = 15.0d0, &
!!$         & h=h, implicit = .true., accuracy_order=1)
!!$    call bdfobj % to_string()
!!$    call bdfobj % solve()
!!$    
!!$    ! find rmse at the last time step
!!$    rmse(2) = 0.0d0
!!$    do j = 1, system % get_num_state_vars()
!!$       do k = 1, bdfobj % num_time_steps
!!$          t = bdfobj % time(k)
!!$          x = system % X(1,j+1)
!!$          error = bdfobj % U (k, 1, j) - exact(t, x)
!!$          rmse(2) = rmse(2) + error**2.0_wp
!!$       end do
!!$    end do
!!$    rmse(2) = sqrt(rmse(2)/(system % get_num_state_vars()*bdfobj % num_time_steps))
!!$
!!$    ! create integrator with supplied mesh resolution
!!$    bdfobj = BDF(system = system, tinit=10.0d0, tfinal = 15.0d0, &
!!$         & h=h, implicit = .true., accuracy_order=2)
!!$    call bdfobj % to_string()
!!$    call bdfobj % solve()
!!$
!!$    ! find rmse at the last time step
!!$    rmse(3) = 0.0d0
!!$    do j = 1, system % get_num_state_vars()
!!$       do k = 1, bdfobj % num_time_steps
!!$          t = bdfobj % time(k)
!!$          x = system % X(1,j+1)
!!$          error = bdfobj % U (k, 1, j) - exact(t, x)
!!$          rmse(3) = rmse(3) + error**2.0_wp
!!$       end do
!!$    end do
!!$    rmse(3) = sqrt(rmse(3)/(system % get_num_state_vars()*bdfobj % num_time_steps))
!!$
!!$    deallocate(system)
!!$
!!$  end subroutine evaluate_temporal_error

  pure real(wp) function exact(t,x) result(val)

    real(wp), intent(in) :: t, x
    real(wp), parameter :: PI = 4.0_wp*atan(1.0_wp)
    real(wp), parameter :: gamma = 0.01_wp

    val = (4.0_wp*pi*gamma*t)**(-0.5_wp)*exp(-((x-t)**2.0_wp)/(4.0_wp*gamma*t))

  end function exact
  
  subroutine evaluate_spatial_error(npts, rmse)

    use dynamic_physics_interface             , only : dynamics
    use unsteady_transport_class              , only : unsteady_transport   

    use abm_integrator_class                  , only : ABM
    use runge_kutta_integrator_class          , only : DIRK
    
    integer         , intent(in)  :: npts
    real(wp)        , intent(out) :: rmse(6)
    type(scalar)    , parameter   :: bounds(2) = [5.0_wp, 45.0_wp]
    class(dynamics) , allocatable :: system
    integer                       :: j, k
    real(wp)                      :: error, t, x
    
    type(ABM)     :: exp_euler, imp_euler, cni
    type(dirk)    :: dirkobj

    ! Case 1
    allocate(system, source = unsteady_transport( &
         & diffusion_coeff = 0.01_WP, &
         & convective_velocity = 1.0_WP, &
         & bounds = bounds, npts=npts, &
         & sparse = .true.))

    ! dirk2
    dirkobj = DIRK(system = system, tinit=10.0d0, tfinal = 11.0d0, &
         & h=1.0d-2, implicit=.true., accuracy_order=2)
    call dirkobj % to_string()
    call dirkobj % solve()
    rmse(1) = 0.0d0
    do j = 1, system % get_num_state_vars()
       do k = 2,2!,1, dirkobj % num_time_steps          
          t = dirkobj % time((dirkobj % num_stages+1)*k - dirkobj % num_stages)
          x = system % X(1,j)
          error = dirkobj % U ((dirkobj % num_stages+1)*k - dirkobj % num_stages, 1, j) - exact(t, x)
          rmse(1) = rmse(1) + error**2.0_wp
       end do
    end do
    rmse(1) = sqrt(rmse(1)/(system % get_num_state_vars()))

    ! dirk 3
    dirkobj = DIRK(system = system, tinit=10.0d0, tfinal = 11.0d0, &
         & h=1.0d-2, implicit=.true., accuracy_order=3)
    call dirkobj % to_string()
    call dirkobj % solve()
    rmse(2) = 0.0d0
    do j = 1, system % get_num_state_vars()
       do k = 2,2!1, dirkobj % num_time_steps          
          t = dirkobj % time((dirkobj % num_stages+1)*k - dirkobj % num_stages)
          x = system % X(1,j)
          error = dirkobj % U ((dirkobj % num_stages+1)*k - dirkobj % num_stages, 1, j) - exact(t, x)
          rmse(2) = rmse(2) + error**2.0_wp
       end do
    end do
    rmse(2) = sqrt(rmse(2)/(system % get_num_state_vars()))

    ! dirk 4
    dirkobj = DIRK(system = system, tinit=10.0d0, tfinal = 11.0d0, &
         & h=1.0d-2, implicit=.true., accuracy_order=4)
    call dirkobj % to_string()
    call dirkobj % solve()
    rmse(3) = 0.0d0
    do j = 1, system % get_num_state_vars()
       do k = 2,2!1, dirkobj % num_time_steps          
          t = dirkobj % time((dirkobj % num_stages+1)*k - dirkobj % num_stages)
          x = system % X(1,j)
          error = dirkobj % U ((dirkobj % num_stages+1)*k - dirkobj % num_stages, 1, j) - exact(t, x)
          rmse(3) = rmse(3) + error**2.0_wp
       end do
    end do
    rmse(3) = sqrt(rmse(3)/(system % get_num_state_vars()))

    ! explicit euler
    exp_euler = ABM(system = system, tinit=10.0d0, tfinal = 11.0d0, &
         & h=1.0d-2, implicit = .false., accuracy_order=1)
    call exp_euler % to_string()
    call exp_euler % solve()
    rmse(5) = 0.0d0
    do j = 1, system % get_num_state_vars()
       do k = 2,2!1, exp_euler % num_time_steps          
          t = exp_euler % time((exp_euler % num_stages+1)*k - exp_euler % num_stages)
          x = system % X(1,j)
          error = exp_euler % U ((exp_euler % num_stages+1)*k - exp_euler % num_stages, 1, j) - exact(t, x)
          rmse(5) = rmse(5) + error**2.0_wp
       end do
    end do
    rmse(5) = sqrt(rmse(5)/(system % get_num_state_vars()))
    
    ! implicit euler
    imp_euler = ABM(system = system, tinit=10.0d0, tfinal = 11.0d0, &
         & h=1.0d-2, implicit = .true., accuracy_order=1)
    call imp_euler % to_string()
    call imp_euler % solve()
    rmse(4) = 0.0d0
    do j = 1, system % get_num_state_vars()
       do k = 2,2!1, imp_euler % num_time_steps          
          t = imp_euler % time((imp_euler % num_stages+1)*k - imp_euler % num_stages)
          x = system % X(1,j)
          error = imp_euler % U ((imp_euler % num_stages+1)*k - imp_euler % num_stages, 1, j) - exact(t, x)
          rmse(4) = rmse(4) + error**2.0_wp
       end do
    end do
    rmse(4) = sqrt(rmse(4)/(system % get_num_state_vars()))

    ! cni
    cni = ABM(system = system, tinit=10.0d0, tfinal = 11.0d0, &
         & h=1.0d-2, implicit = .true., accuracy_order=2)    
    call cni % to_string()
    call cni % solve()
    rmse(6) = 0.0d0
    do j = 1, system % get_num_state_vars()
       do k = 2,2 !1, cni % num_time_steps          
          t = cni % time((cni % num_stages+1)*k - cni % num_stages)
          x = system % X(1,j)
          error = cni % U ((cni % num_stages+1)*k - cni % num_stages, 1, j) - exact(t, x)
          rmse(6) = rmse(6) + error**2.0_wp
       end do
    end do
    rmse(6) = sqrt(rmse(6)/(system % get_num_state_vars()))

    deallocate(system)

  end subroutine evaluate_spatial_error

  subroutine test_integrators(test_system, prefix)

    use abm_integrator_class         , only : ABM
    use runge_kutta_integrator_class , only : DIRK

    class(dynamics) , intent(inout) :: test_system
    character(len=*), intent(in)    :: prefix

    type(ABM)     :: exp_euler, imp_euler, cni
    type(dirk)    :: dirkobj

    exp_euler = ABM(system = test_system, tinit=10.0d0, tfinal = 40.0d0, &
         & h=1.0d-2, implicit=.false., accuracy_order=1)
    call exp_euler % to_string()
    call exp_euler % solve()
    call exp_euler % write_solution(trim(prefix)//"-transport-explicit-euler.dat")

    imp_euler = ABM(system = test_system, tinit=10.0d0, tfinal = 40.0d0, &
         & h=1.0d-2, implicit=.true., accuracy_order=1)
    call imp_euler % to_string()
    call imp_euler % solve()
    call imp_euler % write_solution(trim(prefix)//"-transport-implicit-euler.dat")

    cni = ABM(system = test_system, tinit=10.0d0, tfinal = 40.0d0, &
         & h=1.0d-2, implicit=.true., accuracy_order=2)
    call cni % to_string()
    call cni % solve()
    call cni % write_solution(trim(prefix)//"-transport-cni.dat")

    dirkobj = DIRK(system = test_system, tinit=10.0d0, tfinal = 40.0d0, &
         & h=1.0d-2, implicit=.true., accuracy_order=2)
    call dirkobj % to_string()
    call dirkobj % solve()
    call dirkobj % write_solution(trim(prefix)//"-transport-implicit-dirk2.dat")

    dirkobj = DIRK(system = test_system, tinit=10.0d0, tfinal = 40.0d0, &
         & h=1.0d-2, implicit=.true., accuracy_order=3)
    call dirkobj % to_string()
    call dirkobj % solve()
    call dirkobj % write_solution(trim(prefix)//"-transport-implicit-dirk3.dat")

    dirkobj = DIRK(system = test_system, tinit=10.0d0, tfinal = 40.0d0, &
         & h=1.0d-2, implicit=.true., accuracy_order=4)
    call dirkobj % to_string()
    call dirkobj % solve()
    call dirkobj % write_solution(trim(prefix)//"-transport-implicit-dirk4.dat")

  end subroutine test_integrators

end program test_time_integration

