!=====================================================================!
! Solve one dimensional transport equations using integrators
!=====================================================================!

#include "scalar.fpp"

program test_time_integration

  use constants                 , only : WP
  use dynamic_physics_interface , only : dynamics
  use unsteady_transport_class  , only : unsteady_transport
  
  implicit none

  temporal_convergence : block

    real(wp) :: rmse(3), h
    integer  :: nsteps = 1000
    integer  :: n

    ! grid spacing
    open(12, file='imp-euler-temporal-error.dat')
    write(12, *) "nsteps ", "h ", "rmse1 ", "rmse2 ", "rmse3 "

    do n = 1, 5
       h = 5.0d0/dble(nsteps+1)
       call evaluate_temporal_error(nsteps, rmse)
       write(12, *) nsteps, h, rmse(1), rmse(2), rmse(3)
       nsteps = nsteps*2
    end do

    close(12)

  end block temporal_convergence


  stop

  spatial_convergence : block

    real(wp) :: rmse, h
    integer  :: npts = 500
    integer  :: n

    ! grid spacing
    open(12, file='imp-euler-spatial-error.dat')
    write(12, *) "npts ", "h ", "rmse "

    do n = 1, 5
       h = 40.0d0/dble(npts+1)
       call evaluate_spatial_error(npts, rmse)
       write(12, *) npts, h, rmse
       npts = npts*2
    end do

    close(12)
    
  end block spatial_convergence

!!$
!!$  test_transport: block
!!$
!!$  class(dynamics), allocatable :: system
!!$  type(scalar)   , parameter   :: bounds(2) = [5.0_wp, 45.0_wp]
!!$
!!$    ! case 1
!!$    allocate(system, source = unsteady_transport( &
!!$         & diffusion_coeff = 0.01_WP, &
!!$         & convective_velocity = 1.0_WP, &
!!$         & bounds = bounds, npts=500, &
!!$         & sparse = .true.))
!!$    call test_integrators(system, 'case1')
!!$    deallocate(system)
!!$    
!!$    !case 2
!!$    allocate(system, source = unsteady_transport( &
!!$         & diffusion_coeff = 0.0_WP, &
!!$         & convective_velocity = 1.0_WP, &
!!$         & bounds = bounds, npts=500, &
!!$         & sparse = .true.))
!!$    call test_integrators(system, 'case2')
!!$    deallocate(system)
!!$
!!$  end block test_transport

contains


  subroutine evaluate_temporal_error(nsteps, rmse)

    use dynamic_physics_interface             , only : dynamics
    use backward_differences_integrator_class , only : bdf
    use unsteady_transport_class              , only : unsteady_transport   

    integer         , intent(in)  :: nsteps
    real(wp)        , intent(inout) :: rmse(3)
    type(scalar)    , parameter   :: bounds(2) = [5.0_wp, 45.0_wp]
    class(dynamics) , allocatable :: system
    type(bdf)                     :: bdfobj
    integer                       :: j, k
    real(wp)                      :: error, t, x
    real(wp) :: h

    h = (15.0_wp - 10.0_wp)/dble(nsteps+1)

    ! Case 1
    allocate(system, source = unsteady_transport( &
         & diffusion_coeff = 0.01_WP, &
         & convective_velocity = 1.0_WP, &
         & bounds = bounds, npts = 1000, &
         & sparse = .true.))

    ! create integrator with supplied mesh resolution
    bdfobj = BDF(system = system, tinit=10.0d0, tfinal = 15.0d0, &
         & h=h, implicit = .false., accuracy_order=1)
    call bdfobj % to_string()
    call bdfobj % solve()

    ! find rmse at the last time step
    rmse(1) = 0.0d0   
    loop_vars1 : do j = 1, system % get_num_state_vars()
       loop_time1: do k = 1, bdfobj % num_time_steps
          t = bdfobj % time(k)
          x = system % X(1,j+1)
          error = bdfobj % U (k, 1, j) - exact(t, x)
          rmse(1) = rmse(1) + error**2.0_wp
       end do loop_time1
    end do loop_vars1
    rmse(1) = sqrt(rmse(1)/(system % get_num_state_vars()*bdfobj % num_time_steps))

    ! create integrator with supplied mesh resolution
    bdfobj = BDF(system = system, tinit=10.0d0, tfinal = 15.0d0, &
         & h=h, implicit = .true., accuracy_order=1)
    call bdfobj % to_string()
    call bdfobj % solve()
    
    ! find rmse at the last time step
    rmse(2) = 0.0d0
    loop_vars2 : do j = 1, system % get_num_state_vars()
       loop_time2: do k = 1, bdfobj % num_time_steps
          t = bdfobj % time(k)
          x = system % X(1,j+1)
          error = bdfobj % U (k, 1, j) - exact(t, x)
          rmse(2) = rmse(2) + error**2.0_wp
       end do loop_time2
    end do loop_vars2
    rmse(2) = sqrt(rmse(2)/(system % get_num_state_vars()*bdfobj % num_time_steps))

    ! create integrator with supplied mesh resolution
    bdfobj = BDF(system = system, tinit=10.0d0, tfinal = 15.0d0, &
         & h=h, implicit = .true., accuracy_order=2)
    call bdfobj % to_string()
    call bdfobj % solve()

    ! find rmse at the last time step
    rmse(3) = 0.0d0
    loop_vars3 : do j = 1, system % get_num_state_vars()
       loop_time3: do k = 1, bdfobj % num_time_steps
          t = bdfobj % time(k)
          x = system % X(1,j+1)
          error = bdfobj % U (k, 1, j) - exact(t, x)
          rmse(3) = rmse(3) + error**2.0_wp
       end do loop_time3
    end do loop_vars3
    rmse(3) = sqrt(rmse(3)/(system % get_num_state_vars()*bdfobj % num_time_steps))

    deallocate(system)

  end subroutine evaluate_temporal_error

  pure real(wp) function exact(t,x) result(val)

    real(wp), intent(in) :: t, x
    real(wp), parameter :: PI = 4.0_wp*atan(1.0_wp)
    real(wp), parameter :: gamma = 0.01_wp

    val = (4.0_wp*pi*gamma*t)**(-0.5_wp)*exp(-((x-t)**2.0_wp)/(4.0_wp*gamma*t))

  end function exact

  subroutine evaluate_spatial_error(npts, rmse)

    use dynamic_physics_interface             , only : dynamics
    use backward_differences_integrator_class , only : bdf
    use unsteady_transport_class              , only : unsteady_transport   

    integer         , intent(in)  :: npts
    real(wp)        , intent(out) :: rmse
    type(scalar)    , parameter   :: bounds(2) = [5.0_wp, 45.0_wp]
    class(dynamics) , allocatable :: system
    type(bdf)                     :: bdfobj
    integer                       :: j, k
    real(wp)                      :: error, t, x

    ! Case 1
    allocate(system, source = unsteady_transport( &
         & diffusion_coeff = 0.01_WP, &
         & convective_velocity = 1.0_WP, &
         & bounds = bounds, npts=npts, &
         & sparse = .true.))

    ! create integrator with supplied mesh resolution
    bdfobj = BDF(system = system, tinit=10.0d0, tfinal = 15.0d0, &
         & h=1.0d-4, implicit = .true., accuracy_order=1)
    call bdfobj % to_string()
    call bdfobj % solve()

    ! find rmse at the last time step
    rmse = 0.0d0
    loop_vars : do j = 1, system % get_num_state_vars()
       loop_time: do k = 1, bdfobj % num_time_steps
          t = bdfobj % time(k)
          x = system % X(1,j+1)
          !print *, t, x, exact(t, x), bdfobj % U (k, 1, j)
          error = bdfobj % U (k, 1, j) - exact(t, x)
          rmse = rmse + error**2.0_wp
       end do loop_time
    end do loop_vars
    rmse = sqrt(rmse/(system % get_num_state_vars()*bdfobj % num_time_steps))

    deallocate(system)

  end subroutine evaluate_spatial_error

  subroutine test_integrators(test_system, prefix)

    use abm_integrator_class , only : ABM
    use newmark_integrator_class , only : newmark
    use runge_kutta_integrator_class , only : dirk
    use backward_differences_integrator_class , only : bdf

    class(dynamics) , intent(inout) :: test_system
    character(len=*), intent(in)    :: prefix

    type(ABM)     :: abmobj
    type(newmark) :: nbg
    type(dirk)    :: dirkobj
    type(bdf)     :: bdfobj

!!$
!!$    abmobj = ABM(system = test_system, tinit=10.0d0, tfinal = 20.0d0, &
!!$         & h=1.0d-3, implicit=.false., accuracy_order=1)
!!$    call abmobj % to_string()
!!$    call abmobj % solve()
!!$    call abmobj % write_solution("transport-abm.dat")
!!$
!!$stop
!!$    dirkobj = DIRK(system = test_system, tinit=10.0d0, tfinal = 20.0d0, &
!!$         & h=1.0d-3, implicit=.true., accuracy_order=2)
!!$    call dirkobj % to_string()
!!$    call dirkobj % solve()
!!$    call dirkobj % write_solution("transport-dirk.dat")
    
    bdfobj = BDF(system = test_system, tinit=10.0d0, tfinal = 40.0d0, &
         & h=1.0d-2, implicit = .false., accuracy_order=1)
    call bdfobj % to_string()
    call bdfobj % solve()
    call bdfobj % write_solution(trim(prefix)//"-transport-explicit-euler.dat")   

    bdfobj = BDF(system = test_system, tinit=10.0d0, tfinal = 40.0d0, &
         & h=1.0d-2, implicit = .true., accuracy_order=1)
    call bdfobj % to_string()
    call bdfobj % solve()
    call bdfobj % write_solution(trim(prefix)//"-transport-implicit-euler.dat")   

    bdfobj = BDF(system = test_system, tinit=10.0d0, tfinal = 40.0d0, &
         & h=1.0d-2, implicit = .true., accuracy_order=2)
    call bdfobj % to_string()
    call bdfobj % solve()
    call bdfobj % write_solution(trim(prefix)//"-transport-cni.dat")

  end subroutine test_integrators

end program test_time_integration

