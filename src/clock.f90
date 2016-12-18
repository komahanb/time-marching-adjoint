!=====================================================================!
! Module that manages wall-clock related routines. The clock is
! defined as an object. Type bound procedures defined on the object
! can be used to manage the clock.
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module clock_class

  use constants, only: ZERO

  implicit none

  private
  public :: clock

  !-------------------------------------------------------------------!
  ! Derived type for clock operations. 
  !
  ! The clock can be started and stopped to measure how long different
  ! routines run.
  !
  ! The intrinsic routine system_clock is used to measure time rather
  ! than cpu_time.
  !-------------------------------------------------------------------!

  type clock

     ! object variables
     logical :: running      = .false. ! is clock running?
     integer :: start_counts = 0       ! counts when started
     real(8)    :: elapsed      = 0.      ! total time elapsed in seconds

   contains

     ! type bound procedures
     procedure :: start    
     procedure :: getElapsed
     procedure :: stop
     procedure :: reset

  end type clock

contains

  !-------------------------------------------------------------------!
  ! Start the timer
  !-------------------------------------------------------------------!

  subroutine start(this)
    
    class(clock), intent(inout) :: this

    ! Turn clock on and measure starting time
    this % running = .true.
    call system_clock(this % start_counts)

  end subroutine start

  !-------------------------------------------------------------------!
  ! Start the timer
  !-------------------------------------------------------------------!

  function getElapsed(this) result(elapsed)

    class(clock), intent(in):: this

    integer :: end_counts   ! current number of counts
    integer :: count_rate   ! system-dependent counting rate

    real(8) :: elapsed      ! total elapsed time
    real(8) :: elapsed_time ! elapsed time since last start

    if (this % running) then
       call system_clock(end_counts, count_rate)
       elapsed_time = dble(end_counts - this % start_counts)/dble(count_rate)
       elapsed = this % elapsed + elapsed_time
    else
       elapsed = this % elapsed
    end if

  end function getElapsed

  !-------------------------------------------------------------------!
  ! Stop the clock
  !-------------------------------------------------------------------!  

  subroutine stop(this)

    class(clock), intent(inout) :: this

    ! Check to make sure clock was running
    if (.not. this % running) return

    ! Stop clock and add time
    this % elapsed = this % getElapsed()
    this % running = .false.

  end subroutine stop

  !-------------------------------------------------------------------!
  ! Reset the clock
  !-------------------------------------------------------------------!

  subroutine reset(this)

    class(clock), intent(inout) :: this

    this % running      = .false.
    this % start_counts = 0
    this % elapsed      = ZERO

  end subroutine reset

end module clock_class
