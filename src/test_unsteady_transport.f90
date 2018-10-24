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

  convergence : block

    class(dynamics), allocatable :: system
    type(scalar)   , parameter   :: bounds(2) = [5.0_wp, 45.0_wp]
    real(wp) :: rmse(6), time(6), h
    integer  :: npts
    integer  :: n

    !-----------------------------------------------------------------!
    ! Run for convection - diffusion physics
    !-----------------------------------------------------------------!

    ! File handlers
    open(13, file='timing-case1.dat')
    write(13, *) "npts ", "h ", &
         & "dirk2 ", "dirk3 ", "dirk4 ", &
         & "imp-euler ", "exp-euler ", "cni "
    
    open(12, file='spatial-error-case1.dat')
    write(12, *) "npts ", "h ", &
         & "dirk2 ", "dirk3 ", "dirk4 ", &
         & "imp-euler ", "exp-euler ","cni "

    ! Reference spacing
    h   = 0.1d0
    npts = 400

    ! Solve for each spacing and find rmse and time
    do n = 1, 5

       allocate(system, source = unsteady_transport( &
            & diffusion_coeff = 0.01_WP, &
            & convective_velocity = 1.0_WP, &
            & bounds = bounds, npts = npts, &
            & sparse = .true.))
       
       call evaluate_spatial_error(system, exact_conv_diff, h, rmse, time)
       
       write(12, *) npts, h, rmse(1), rmse(2), rmse(3), rmse(4), rmse(5), rmse(6)
       write(13, *) npts, h, time(1), time(2), time(3), time(4), time(5), time(6)

       ! Refine spacing
       h = h/2.0d0
       npts = npts*2

       deallocate(system)

    end do

    close(12)
    close(13)
   
    !-----------------------------------------------------------------!
    ! Run for diffusion physics
    !-----------------------------------------------------------------!
    
    ! File handlers
    open(13, file='timing-case2.dat')
    write(13, *) "npts ", "h ", &
         & "dirk2 ", "dirk3 ", "dirk4 ", &
         & "imp-euler ", "exp-euler ", "cni "

    open(12, file='spatial-error-case2.dat')
    write(12, *) "npts ", "h ", &
         & "dirk2 ", "dirk3 ", "dirk4 ", &
         & "imp-euler ", "exp-euler ","cni "

    ! Reference spacing
    h = 0.1d0
    npts = 400

    ! Solve for each spacing and find rmse and time
    do n = 1, 5

       allocate(system, source = unsteady_transport( &
            & diffusion_coeff = 0.00_WP, &
            & convective_velocity = 1.0_WP, &
            & bounds = bounds, npts = npts, &
            & sparse = .true.))

       call evaluate_spatial_error(system, exact_diff, h, rmse, time)

       write(12, *) npts, h, rmse(1), rmse(2), rmse(3), rmse(4), rmse(5), rmse(6)
       write(13, *) npts, h, time(1), time(2), time(3), time(4), time(5), time(6)

       ! Refine spacing
       h = h/2.0d0
       npts = npts*2

       deallocate(system)

    end do

    close(12)
    close(13)

  end block convergence

contains

  subroutine test_integrators(test_system, prefix)

    use abm_integrator_class         , only : ABM
    use runge_kutta_integrator_class , only : DIRK
    use backward_differences_integrator_class, only : BDF

    class(dynamics) , intent(inout) :: test_system
    character(len=*), intent(in)    :: prefix

    type(BDF)     :: exp_euler
    type(ABM)     :: imp_euler, cni
    type(dirk)    :: dirkobj
  
    exp_euler = BDF(system = test_system, tinit=10.0d0, tfinal = 40.0d0, &
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


  ! Exact solution for convection-diffusion  
  pure real(wp) function exact_conv_diff(t,x) result(val)

    real(wp), intent(in) :: t, x
    real(wp), parameter :: PI = 4.0_wp*atan(1.0_wp)
    real(wp), parameter :: gamma = 0.01_wp

    val = (4.0_wp*pi*gamma*t)**(-0.5_wp)*exp(-((x-t)**2.0_wp)/(4.0_wp*gamma*t))   
    
  end function exact_conv_diff

  ! Exact solution for diffusion
  pure real(wp) function exact_diff(t,x) result(val)

    real(wp), intent(in) :: t, x
    real(wp), parameter :: PI = 4.0_wp*atan(1.0_wp)
    real(wp), parameter :: gamma = 0.01_wp

    val = (0.4_wp*PI)**(-0.5_wp)*exp(-2.5_wp*(x-t)*(x-t))
    
  end function exact_diff

  subroutine evaluate_spatial_error(system, exact, h, rmse, time)

    use dynamic_physics_interface             , only : dynamics
    use unsteady_transport_class              , only : unsteady_transport    
    use backward_differences_integrator_class , only : bdf
    use abm_integrator_class                  , only : ABM
    use runge_kutta_integrator_class          , only : DIRK
    
    ! Arguments
    class(dynamics), intent(in) :: system   
    interface
       pure real(8) function exact(t, x)
         real(8), intent(in) ::  t, x
       end function exact
    end interface
    real(wp), intent(in) :: h
    real(wp) , intent(out) :: rmse(6), time(6)

    ! Integrators used
    type(ABM)     :: imp_euler, cni
    type(BDF)     :: exp_euler
    type(dirk)    :: dirkobj

    ! Locals
    real(wp) :: tic, toc
    integer :: j, k
    real(wp) :: error, t, x
    
    ! dirk2
    call cpu_time(tic)           
    dirkobj = DIRK(system = system, tinit=10.0d0, tfinal = 40.0d0, &
         & h=h, implicit=.true., accuracy_order=2)
    call dirkobj % to_string()
    call dirkobj % solve()
    call cpu_time(toc)
    time(1) = toc - tic
    rmse(1) = 0.0d0
    do j = 1, system % get_num_state_vars()
       do k = 1, dirkobj % num_time_steps          
          t = dirkobj % time((dirkobj % num_stages+1)*k - dirkobj % num_stages)
          x = system % X(1,j)
          error = dirkobj % U ((dirkobj % num_stages+1)*k - dirkobj % num_stages, 1, j) - exact(t, x)
          rmse(1) = rmse(1) + error**2.0_wp
       end do
    end do
    rmse(1) = sqrt(rmse(1)/(dirkobj % num_time_steps * system % get_num_state_vars()))

    ! dirk 3
    call cpu_time(tic)
    dirkobj = DIRK(system = system, tinit=10.0d0, tfinal = 40.0d0, &
         & h=h, implicit=.true., accuracy_order=3)
    call dirkobj % to_string()
    call dirkobj % solve()
    call cpu_time(toc)
    time(2) = toc - tic
    rmse(2) = 0.0d0
    do j = 1, system % get_num_state_vars()
       do k = 1, dirkobj % num_time_steps          
          t = dirkobj % time((dirkobj % num_stages+1)*k - dirkobj % num_stages)
          x = system % X(1,j)
          error = dirkobj % U ((dirkobj % num_stages+1)*k - dirkobj % num_stages, 1, j) - exact(t, x)
          rmse(2) = rmse(2) + error**2.0_wp
       end do
    end do
    rmse(2) = sqrt(rmse(2)/(dirkobj % num_time_steps * system % get_num_state_vars()))

    ! dirk 4
    call cpu_time(tic)
    dirkobj = DIRK(system = system, tinit=10.0d0, tfinal = 40.0d0, &
         & h=h, implicit=.true., accuracy_order=4)
    call dirkobj % to_string()
    call dirkobj % solve()
    call cpu_time(toc)
    time(3) = toc - tic
    rmse(3) = 0.0d0
    do j = 1, system % get_num_state_vars()
       do k = 1, dirkobj % num_time_steps          
          t = dirkobj % time((dirkobj % num_stages+1)*k - dirkobj % num_stages)
          x = system % X(1,j)
          error = dirkobj % U ((dirkobj % num_stages+1)*k - dirkobj % num_stages, 1, j) - exact(t, x)
          rmse(3) = rmse(3) + error**2.0_wp
       end do
    end do
    rmse(3) = sqrt(rmse(3)/(dirkobj % num_time_steps * system % get_num_state_vars()))
    
    ! implicit euler
    call cpu_time(tic)
    imp_euler = ABM(system = system, tinit=10.0d0, tfinal = 40.0d0, &
         & h=h, implicit = .true., accuracy_order=1)
    call imp_euler % to_string()
    call imp_euler % solve()
    call cpu_time(toc)
    time(4) = toc - tic
    rmse(4) = 0.0d0
    do j = 1, system % get_num_state_vars()
       do k = 1, imp_euler % num_time_steps          
          t = imp_euler % time((imp_euler % num_stages+1)*k - imp_euler % num_stages)
          x = system % X(1,j)
          error = imp_euler % U ((imp_euler % num_stages+1)*k - imp_euler % num_stages, 1, j) - exact(t, x)
          rmse(4) = rmse(4) + error**2.0_wp
       end do
    end do
    rmse(4) = sqrt(rmse(4)/(imp_euler % num_time_steps * system % get_num_state_vars()))

    ! explicit euler
    call cpu_time(tic)
    exp_euler = BDF(system = system, tinit=10.0d0, tfinal = 40.0d0, &
         & h=h, implicit = .false., accuracy_order=1)
    call exp_euler % to_string()
    call exp_euler % solve()
    call cpu_time(toc)
    time(5) = toc - tic
    rmse(5) = 0.0d0
    do j = 1, system % get_num_state_vars()
       do k = 1, exp_euler % num_time_steps          
          t = exp_euler % time((exp_euler % num_stages+1)*k - exp_euler % num_stages)
          x = system % X(1,j)
          error = exp_euler % U ((exp_euler % num_stages+1)*k - exp_euler % num_stages, 1, j) - exact(t, x)
          rmse(5) = rmse(5) + error**2.0_wp
       end do
    end do
    rmse(5) = sqrt(rmse(5)/(exp_euler % num_time_steps*system % get_num_state_vars()))
    
    ! cni
    call cpu_time(tic)
    cni = ABM(system = system, tinit=10.0d0, tfinal = 40.0d0, &
         & h=h, implicit = .true., accuracy_order=2)    
    call cni % to_string()
    call cni % solve()
    call cpu_time(toc)
    time(6) = toc - tic
    rmse(6) = 0.0d0
    do j = 1, system % get_num_state_vars()
       do k = 1, cni % num_time_steps          
          t = cni % time((cni % num_stages+1)*k - cni % num_stages)
          x = system % X(1,j)
          error = cni % U ((cni % num_stages+1)*k - cni % num_stages, 1, j) - exact(t, x)
          rmse(6) = rmse(6) + error**2.0_wp
       end do
    end do
    rmse(6) = sqrt(rmse(6)/(cni % num_time_steps*system % get_num_state_vars()))

  end subroutine evaluate_spatial_error

end program test_time_integration

