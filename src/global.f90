#include "scalar.fpp"

!---------------------------------------------------------------------!
!  Module containing all global variables used in the simulation
!---------------------------------------------------------------------!

module global

  use header, only: heat_obj
  use timing, only: timer

  implicit none
  save

  ! Main object
  type(heat_obj) :: heat

  ! Timing objects
  type(timer) :: time_total ! timer for whole calculation
  type(timer) :: time_mat   ! timer for mat building
  type(timer) :: time_res   ! timer for res building
  type(timer) :: time_lin   ! timer for linear solves

  ! PETSc error code
  type(integer) :: ierr

  ! MPI parameters
  type(logical) :: master = .false. ! am I master?
  type(integer) :: rank             ! rank of processor
  type(integer) :: nprocs           ! number of processors

  ! Solver type
  character(len=25) :: solver_type

end module global
