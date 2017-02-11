!=====================================================================!
! Module to handle all MPI tasks
!
! Author: Komahan Boopathy (komahan@gatech.edu)
!=====================================================================!

module mpi_wrapper
  
  use constants, only: LOG_UNIT
  use variables, only: idproc, numproc, master, mpierr

#if defined USE_COMPLEX
  use mpi
#endif
  implicit none
  
  private
  public :: MPI_START_ALL, MPI_STOP_ALL
  
contains

  !-------------------------------------------------------------------!
  ! Perform MPI initialization, identify master and say hello!
  !-------------------------------------------------------------------!

  subroutine MPI_START_ALL()

#if defined MPI
    ! Initialize MPI
    call MPI_Init(mpierr)

    if (mpierr /= MPI_SUCCESS) then
       stop 'MPI Initializatin error'
    endif

    ! Get the ID and number of processors
    call MPI_Comm_rank(MPI_COMM_WORLD, idproc , mpierr)
    call MPI_Comm_size(MPI_COMM_WORLD, numproc, mpierr)

    ! find master
    if (idproc .eq. 0) master = .true.

    ! Say hello
    !if(master) then
    !   write(LOG_UNIT,'(2x,a,i3)') '>> Number of Processors = ', numproc
    !end if
    ! write(LOG_UNIT,'(2x,a,i3)') '>> Processor ', idproc + 1, 'of', numproc,'...[OK]'
#endif

  end subroutine MPI_START_ALL

  !-------------------------------------------------------------------!
  ! Stop all the MPI processors                                       !
  !-------------------------------------------------------------------!

  subroutine MPI_STOP_ALL()
#if defined MPI
    ! Make sure all procs are here
    call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    if (mpierr /= MPI_SUCCESS) then
       stop 'MPI Barrier error'
    endif

    ! Finalize MPI
    call MPI_Finalize(mpierr)
    if (mpierr /= MPI_SUCCESS) then
       stop 'MPI Finalize error'
    endif
#endif
  end subroutine MPI_STOP_ALL

end module mpi_wrapper
