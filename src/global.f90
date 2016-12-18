#include "scalar.fpp"

!---------------------------------------------------------------------!
!  Module containing all global variables used in the simulation
!---------------------------------------------------------------------!

module global

  use diffusion_class, only: diffusion
  use clock_class, only: clock

  implicit none

  ! Make all variables static (even if the module is descoped the
  ! values retain their values)
  save

  ! Main object
  type(diffusion) :: diff_obj

  ! Timing objects
  type(clock) :: clock_total ! clock for whole calculation
  type(clock) :: clock_mat   ! clock for mat building
  type(clock) :: clock_res   ! clock for res building
  type(clock) :: clock_lin   ! clock for linear solves

  ! PETSc error code
  type(integer) :: ierr

  ! MPI parameters
  type(logical) :: master = .false. ! am I master?
  type(integer) :: rank             ! rank of processor
  type(integer) :: nprocs           ! number of processors

  ! Solver type
  character(len=25) :: solver_type

end module global
