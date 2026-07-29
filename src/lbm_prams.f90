#ifndef LBM_NX
#define LBM_NX 128
#endif
#ifndef LBM_NY
#define LBM_NY 128
#endif

module lbm_params
  use iso_fortran_env, only: dp => real64
  implicit none

  ! Shared D2Q9 constants in lattice units.
  integer, parameter :: q  = 9
  integer, parameter :: nx = LBM_NX
  integer, parameter :: ny = LBM_NY

  real(dp), parameter :: default_omega = 1.0_dp / 0.7_dp
  real(dp), parameter :: cs2 = 1.0_dp / 3.0_dp

  ! Discrete velocities, opposite directions, and quadrature weights.
  integer, parameter :: cx(q) = [0, 1, 0,-1, 0, 1,-1,-1, 1]
  integer, parameter :: cy(q) = [0, 0, 1, 0,-1, 1, 1,-1,-1]

  integer, parameter :: opp(q) = [1,4,5,2,3,8,9,6,7]

  real(dp), parameter :: w(q) = [ &
       4.0_dp/9.0_dp, &
       1.0_dp/9.0_dp, 1.0_dp/9.0_dp, 1.0_dp/9.0_dp, 1.0_dp/9.0_dp, &
       1.0_dp/36.0_dp, 1.0_dp/36.0_dp, 1.0_dp/36.0_dp, 1.0_dp/36.0_dp ]


  ! Boundary identifiers and the runtime configuration of global edges.
  integer, parameter :: BC_PERIODIC    = 1
  integer, parameter :: BC_WALL        = 2
  integer, parameter :: BC_MOVING_WALL = 3

  integer :: bc_left, bc_right, bc_bottom, bc_top

end module lbm_params
