!=====================================================================!
! Program that solves the heat equation
!=====================================================================!

program main
  
  use constants
  use variables

  implicit none

  type(diffusion) :: hello

  call initialize()
   
  call execute()
 
  call finalize()

contains

  !-------------------------------------------------------------------!
  ! Initialize the simulation
  !-------------------------------------------------------------------!
  
  subroutine initialize()

    use dimmpi
  
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
  ! Execute the simulation
  !-------------------------------------------------------------------!
  
  subroutine execute()

    write(*, *) 'Executing'

  end subroutine execute

  !-------------------------------------------------------------------!
  ! Finalize the simulation
  !-------------------------------------------------------------------!

  subroutine finalize()
    
    use dimmpi

    ! stop timer
    call clock_tot % stop()

    ! print output
    if (master) then
       write(LOG_UNIT, '(/,/,A)') 'RESULTS'
       write(LOG_UNIT, '(A)')     '**********************'
       write(LOG_UNIT, '(/,"Total Execution time (s): ",F0.4)') clock_tot % elapsed
    end if

    call MPI_STOP_ALL

  end subroutine finalize

end program main
