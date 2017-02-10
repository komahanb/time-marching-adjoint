#include "scalar.fpp"

!=====================================================================!
! Program that solves the diffusion equation
!=====================================================================!

program main
  
  use constants
  use variables

  implicit none

  call initialize()
   
  call execute()
 
  call finalize()

contains

  !-------------------------------------------------------------------!
  ! Execute the simulation on one dimensional grid
  !-------------------------------------------------------------------!
  
  subroutine execute()
    
    use mpi_wrapper

    type(integer), parameter :: NDIM    = 1  ! physical dimension
    type(integer), parameter :: NPOINTS = 99  ! number of grid points
    type(integer), parameter :: NSTEPS  = 99  ! number of time steps
    
    ! state vector
    type(scalar), allocatable, dimension(:,:) :: state

    ! diffusion operator matrix
    type(scalar), allocatable, dimension(:,:) :: D  
    
    ! grid
    type(scalar), allocatable, dimension(:,:) :: xpts
    
    ! diffusion parameter
    type(scalar), parameter :: alpha = 1.0d-3
    
    ! generation of grid and time steps
    type(scalar) :: dx, dt

    ! Loop over
    type(integer) :: i, j, k
    type(integer) :: n, p
    type(integer) :: irow, icol

    ! spatial and temporal bounds
    type(scalar), parameter :: xfinal = 1.0_WP
    type(scalar), parameter :: tfinal = 1.0_WP
    
    type(scalar) :: dval, odval

    ! Lapack variables
    integer, allocatable, dimension(:)        :: ipiv
    integer                                   :: info, size

    ! Temporary lapack vector for RHS
    type(scalar), allocatable, dimension(:) :: rhs
    
    allocate(state(NPOINTS, NSTEPS))
    allocate(rhs(NPOINTS))
    allocate(xpts(NDIM,NPOINTS))    
    allocate(D(NPOINTS, NPOINTS))
    
    dx = xfinal/dble(NPOINTS-1)
    dt = tfinal/dble(NSTEPS-1)
    
    ! generate the spatial grid
    do i = 1, NPOINTS
       xpts(:,i) = dble(i-1)*dx
    end do

    ! Initial condition
    do p = 1, NPOINTS       
       state(p,1) = sin(10.0_WP*PI*xpts(1,p)) ! note that
    end do
   
    ! Boundary conditions

    ! left end
    state(1, :) = 0.15_WP

    ! right end
    state(NPOINTS, :) = 0.15_WP

    odval = alpha*dt/(dx*dx)
    dval  = 1.0_WP + 2.0_WP*odval
    
    ! Assemble diffusion matrix using central difference
    D = 0.0_WP
    do irow = 1, NPOINTS

       ! leading diagonal
       icol = irow
       D(irow,icol) = dval

       ! sub diagonal
       icol = irow - 1
       if ( icol .ge. 1) then
          D(irow, icol) = -odval
       end if

       ! super diagonal
       icol = irow + 1
       if ( icol .le. NPOINTS) then
          D(irow, icol) = -odval
       end if

    end do

    ! apply BC to the matrix ( make the diagonal as identity, and off
    ! diagonal as zero)
    ! off diagonal as zero
    D(1,:) = 0.0_WP
    D(NPOINTS,:) = 0.0_WP

    D(1,1) = 1.0_WP
    D(NPOINTS,NPOINTS) = 1.0_WP

    !D =  transpose(D)
  
!!$    do i = 1, NPOINTS
!!$       print *, i, D(i,:)
!!$    end do

    size = NPOINTS
    if ( .not. allocated(ipiv)   ) allocate( ipiv(size) )

    do n = 2, NSTEPS
       
       ! Previous solution is the RHS
       rhs = state(:, n-1)                    
       
#if defined COMPLEX
       call ZGESV(size, 1, D, size, IPIV, rhs, size, INFO)
#else
       call DGESV(size, 1, D, size, IPIV, rhs, size, INFO)
#endif
       if (INFO .ne. 0) then
          print*, "LAPACK ERROR:", info
          call MPI_STOP_ALL
       end if
       
       ! Result is the next solution
       state(:, n)  = rhs

       write(*,'(F15.3,E15.3)') (n-1)*dt, norm2(state(:,n-1) - state(:,n))/sqrt(dble(NPOINTS))

    end do

    ! write the solution
    open (unit = 2, file = "solution.dat")
    write(2,*) 'variables= "x","t","u"'
    write(2,*) 'zone i=', NPOINTS,' j=',NSTEPS,' datapacking=point'
    do i = 1, NSTEPS
       do j = 1, NPOINTS
          write (2,*) dble(j-1)*dx, dble(i-1)*dt, state(j,i)
       end do
    end do
    close(2)

    if (allocated(ipiv)) deallocate(ipiv)

    deallocate(rhs)
    deallocate(state)
    deallocate(xpts)
    deallocate(D)

  end subroutine execute

  !-------------------------------------------------------------------!
  ! Initialize the simulation
  !-------------------------------------------------------------------!
  
  subroutine initialize()

    use mpi_wrapper
  
    integer :: i                    ! loop counter
    integer :: argc                 ! number of command line arguments
    character(len=25) :: argv       ! a command line argument
    character(len=10) :: today_date
    character(len=8)  :: today_time

    ! INITIALIZE MPI   

    call MPI_START_ALL()

    ! PRINT DETAILS

    if (master) then
       
       write(LOG_UNIT,*)
       write(LOG_UNIT,*) '======================================================='
       write(LOG_UNIT,*) '                    DIFFUSION                          '
       write(LOG_UNIT,*) '======================================================='
       write(LOG_UNIT,*)

       ! author
       write(LOG_UNIT, *) &
            '           Author:  Komahan Boopathy'

       ! email
       write(LOG_UNIT, *) &
            '            Email:  komahan@gatech.edu'
       
       ! institute
       write(LOG_UNIT, *) &
            '     Developed At:  Georgia Institute of Technology'

       ! version
       write(LOG_UNIT, '(6X,"Version:",7X,I1,".",I1,".",I1)') &
            VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE
       
       ! date and time of execution
       call clock_tot % getDateTime(today_date, today_time)
       write(LOG_UNIT, '(6X,"Date/Time:",5X,A,1X,A)') &
            trim(today_date), trim(today_time)
       
       ! number of processors
       write(LOG_UNIT, '(1X,A,I0)') '     MPI Processes: ',numproc
       
    end if

    ! PROCESS COMMAND LINE ARGS

    ! get number of command line arguments
    argc = COMMAND_ARGUMENT_COUNT()

    ! begin loop around command line args
    do i = 1,argc

      ! get that argument
      call GET_COMMAND_ARGUMENT(i,argv)

      ! begin case structure
      select case(trim(argv))

        ! solver
        case('--solver_type')

          ! get next argument
          call GET_COMMAND_ARGUMENT(i+1,argv)

          ! set global var
          solver_type = trim(argv)

        ! do nothing here
        case DEFAULT

      end select

    end do
    
    ! Start the total timer clock
    call clock_tot % start()

  end subroutine initialize

  !-------------------------------------------------------------------!
  ! Finalize the simulation
  !-------------------------------------------------------------------!

  subroutine finalize()
    
    use mpi_wrapper

    ! stop timer
    call clock_tot % stop()

    ! print output
    if (master) then
       write(LOG_UNIT, '(/,A)') '======================================================='
       write(LOG_UNIT, '(A)')   '                    RESULTS                            '
       write(LOG_UNIT, '(A,/)') '======================================================='
       write(LOG_UNIT, '("Total Execution time (s): ",F0.4)') clock_tot % elapsed
    end if

    call MPI_STOP_ALL

  end subroutine finalize

end program main
