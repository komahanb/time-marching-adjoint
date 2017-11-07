#include "scalar.fpp"

!---------------------------------------------------------------------!
!  Module containing all global variables used in the simulation
!---------------------------------------------------------------------!

module variables

  !use diffusion_class, only: diffusion
  use clock_class, only: clock

  implicit none

  ! Make all variables static (even if the module is descoped the
  ! values retain their values)
  save

  ! Main object
  !type(diffusion) :: physics

  ! Timing objects
  type(clock) :: clock_tot   ! clock for whole calculation
  type(clock) :: clock_mat   ! clock for mat building
  type(clock) :: clock_res   ! clock for res building
  type(clock) :: clock_lin   ! clock for linear solves

  ! PETSc error code
  type(integer) :: ierr

  ! MPI parameters
  type(logical) :: master = .false. ! am I master?
  type(integer) :: idproc           ! rank of processor
  type(integer) :: numproc          ! number of processors
  type(integer) :: mpierr           ! error code for mpi calls

  ! Solver type
  character(len=25) :: solver_type

end module variables
