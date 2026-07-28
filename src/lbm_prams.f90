#ifndef LBM_NX
#define LBM_NX 128
#endif
#ifndef LBM_NY
#define LBM_NY 128
#endif

module lbm_params
  use iso_fortran_env, only: dp => real64
  implicit none

  ! Shared numerical and physical constants for the D2Q9 lattice-Boltzmann
  ! solver.  Lattice units are used throughout: dx = dt = 1.

  ! ---------- Global lattice ----------
  ! q is the number of populations stored at every lattice node: one rest
  ! population, four axis-aligned populations, and four diagonal populations.
  integer, parameter :: q  = 9
  ! Override at compile time with -DLBM_NX=<n> -DLBM_NY=<n>.
  integer, parameter :: nx = LBM_NX
  integer, parameter :: ny = LBM_NY

  ! ---------- Physics ----------
  ! For the BGK collision operator, nu = cs2*(1/omega - 1/2).
  real(dp), parameter :: default_omega = 1.0_dp / 0.7_dp
  ! Squared isothermal lattice speed of sound for D2Q9.
  real(dp), parameter :: cs2 = 1.0_dp / 3.0_dp

  ! ---------- D2Q9 velocities ----------
  ! Direction numbering used throughout the solver:
  !
  !       7  3  6
  !        \ | /
  !     4 -- 1 -- 2       1 is the rest population
  !        / | \
  !       8  5  9
  !
  ! cx/cy are integer displacements because one population moves exactly one
  ! lattice cell per time step.
  integer, parameter :: cx(q) = [0, 1, 0,-1, 0, 1,-1,-1, 1]
  integer, parameter :: cy(q) = [0, 0, 1, 0,-1, 1, 1,-1,-1]

  ! opp maps a population to the direction with the reverse velocity.  This is
  ! the conceptual mapping used by bounce-back boundary conditions.
  integer, parameter :: opp(q) = [1,4,5,2,3,8,9,6,7]

  ! Standard D2Q9 quadrature weights: rest, axial, then diagonal.
  real(dp), parameter :: w(q) = [ &
       4.0_dp/9.0_dp, &
       1.0_dp/9.0_dp, 1.0_dp/9.0_dp, 1.0_dp/9.0_dp, 1.0_dp/9.0_dp, &
       1.0_dp/36.0_dp, 1.0_dp/36.0_dp, 1.0_dp/36.0_dp, 1.0_dp/36.0_dp ]


  ! ---------- Boundary-condition identifiers ----------
  ! Periodic boundaries are implemented by making the image at the opposite
  ! edge a neighbour during halo exchange.  Stationary and moving walls are
  ! physical boundaries and are imposed by bounce-back after streaming.
  integer, parameter :: BC_PERIODIC    = 1
  integer, parameter :: BC_WALL        = 2
  integer, parameter :: BC_MOVING_WALL = 3

  ! Runtime configuration of the four global edges.  Every coarray image gets
  ! the same values before the decomposition and neighbour map are constructed.
  integer :: bc_left, bc_right, bc_bottom, bc_top

end module lbm_params
