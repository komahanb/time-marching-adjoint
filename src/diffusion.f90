#include "scalar.fpp"

!=====================================================================!
! Module for modeling the physics as diffusion
!=====================================================================!

module diffusion_class
  
  use physics_interface, only: physics
  
  ! type, extends(physics) :: diffusion

  type :: diffusion

     ! diffusion coefficient
     type(scalar) :: diffcof

  end type diffusion

  contains




end module diffusion_class
