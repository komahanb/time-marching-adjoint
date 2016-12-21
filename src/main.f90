!=====================================================================!
! Program that solves the heat equation
!=====================================================================!

program heat

  use constants
  use variables

  ! petsc modules
  use petscsys
  use petscmat
  use petscvec

  implicit none

  print *, "Hello world"
 
  ! initialize
  call initialize()
  
  ! begin total timer
  call clock_tot % start()
  
  ! execute cmfd
  call execute()
  
  ! stop timer
  call clock_tot % stop()

  ! finalize run
  call finalize()

contains

  !-------------------------------------------------------------------!
  ! Initialize the simulation
  !-------------------------------------------------------------------!
  
  subroutine initialize()

    integer :: i                    ! loop counter
    integer :: argc                 ! number of command line arguments
    character(len=25) :: argv       ! a command line argument
    character(len=10) :: today_date
    character(len=8)  :: today_time

     ! INITIALIZE MPI
    
    ! initialize petsc
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    
    ! get mpi info
    call MPI_COMM_RANK(PETSC_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(PETSC_COMM_WORLD, nprocs, ierr)
    
    ! find master
    if (rank == 0) master = .true.

    !=========================================
    ! WRITE HEADING INFO

    if (master) then

      ! write header
      write(*, FMT='(/11(A/))') &
    & '    ,o888888o.           ,8.       ,8.          8 8888888888   8 888888888o.      ', &
    & '   8888     `88.        ,888.     ,888.         8 8888         8 8888    `^888.   ', &
    & ',8 8888       `8.      .`8888.   .`8888.        8 8888         8 8888        `88. ', &
    & '88 8888               ,8.`8888. ,8.`8888.       8 8888         8 8888         `88 ', &
    & '88 8888              ,8"8.`8888,8^8.`8888.      8 888888888888 8 8888          88 ', &
    & '88 8888             ,8" `8.`8888" `8.`8888.     8 8888         8 8888          88 ', &
    & '88 8888            ,8"   `8.`88"   `8.`8888.    8 8888         8 8888         ,88 ', &
    & '`8 8888       .8" ,8"     `8.`"     `8.`8888.   8 8888         8 8888        ,88" ', &
    & '   8888     ,88" ,8"       `8        `8.`8888.  8 8888         8 8888    ,o88P"   ', &
    & '    `8888888P"  ,8"         `         `8.`8888. 8 8888         8 888888888P"      ', &
    & '__________________________________________________________________________________'

      ! Write version information
      write(*, FMT=*) &
           '     Developed At:  Georgia Institute of Technology'
      write(*, FMT='(6X,"Version:",7X,I1,".",I1,".",I1)') &
           VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE
      
      ! Write the date and time
      call get_today(today_date, today_time)
      write(*, FMT='(6X,"Date/Time:",5X,A,1X,A)') &
           trim(today_date), trim(today_time)
      
      ! Write information on number of processors
      write(*, FMT='(1X,A,I0)') '     MPI Processes: ',nprocs

    end if

    !=======================================
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

  end subroutine initialize

  !-------------------------------------------------------------------!
  ! Execute the simulation
  !-------------------------------------------------------------------!
  
  subroutine execute()

  end subroutine execute

  !-------------------------------------------------------------------!
  ! Finalize the simulation
  !-------------------------------------------------------------------!

  subroutine finalize()
    
    ! finalize petsc
    call PetscFinalize(ierr)
    
    ! print output
    if (master) then
       write(*, FMT='(/,/,A)') 'RESULTS'
       write(*, FMT='(A)')     '**********************'
       write(*, FMT='(/,"Total Execution time (s): ",F0.4)') clock_tot % elapsed
       ! write(*, FMT='("K-effective: ",F0.8,/)') cmfd%keff
    end if

  end subroutine finalize

!---------------------------------------------------------------------!
! GET_TODAY determines the date and time at which the program began
! execution and returns it in a readable format
! ---------------------------------------------------------------------!

  subroutine get_today(today_date, today_time)

    character(10), intent(out) :: today_date
    character(8),  intent(out) :: today_time
    
    integer       :: val(8)
    character(8)  :: date_
    character(10) :: time_
    character(5)  :: zone
    
    call date_and_time(date_, time_, zone, val)
    ! val(1) = year (YYYY)
    ! val(2) = month (MM)
    ! val(3) = day (DD)
    ! val(4) = timezone
    ! val(5) = hours (HH)
    ! val(6) = minutes (MM)
    ! val(7) = seconds (SS)
    ! val(8) = milliseconds

    if (val(2) < 10) then
       if (val(3) < 10) then
          today_date = date_(6:6) // "/" // date_(8:8) // "/" // date_(1:4)
       else
          today_date = date_(6:6) // "/" // date_(7:8) // "/" // date_(1:4)
       end if
    else
       if (val(3) < 10) then
          today_date = date_(5:6) // "/" // date_(8:8) // "/" // date_(1:4)
       else
          today_date = date_(5:6) // "/" // date_(7:8) // "/" // date_(1:4)
       end if
    end if
    today_time = time_(1:2) // ":" // time_(3:4) // ":" // time_(5:6)

  end subroutine get_today

end program heat
