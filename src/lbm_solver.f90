module lbm_solver
  use lbm_params
  use domain_decomposition
  implicit none
  private

  public :: configure_solver
  public :: initialize_uniform, initialize_shear_wave, initialize_density_wave
  public :: initialize_equilibrium, macro_val, collide, stream
  public :: apply_boundaries, velocity_residual

  ! Single-relaxation-time (BGK) collision frequency.  It is configured once
  ! from the command line and controls the lattice viscosity.
  real(dp) :: relaxation_omega

contains

  subroutine configure_solver(omega)
    real(dp), intent(in) :: omega

    relaxation_omega = omega
  end subroutine configure_solver
  subroutine initialize_uniform(rho_i, ux_i, uy_i)
    real(dp), intent(out) :: rho_i(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(out) :: ux_i(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(out) :: uy_i(0:nx_alloc+1,0:ny_alloc+1)

    ! Initialize the complete allocation, including ghosts and unused padding,
    ! to benign values.  Only 1:nx_loc,1:ny_loc represents owned fluid cells.
    rho_i = 1.0_dp
    ux_i = 0.0_dp
    uy_i = 0.0_dp
  end subroutine initialize_uniform


  subroutine initialize_shear_wave(rho_i, ux_i, uy_i, amplitude)
    real(dp), intent(out) :: rho_i(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(out) :: ux_i(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(out) :: uy_i(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(in) :: amplitude
    integer :: ii, jj, global_y
    real(dp) :: phase
    real(dp), parameter :: pi = acos(-1.0_dp)

    ! Use the global y coordinate so the sinusoid remains continuous across
    ! image boundaries and is independent of the chosen decomposition.
    call initialize_uniform(rho_i,ux_i,uy_i)
    do jj = 1, ny_loc
      global_y = y_offset + jj - 1
      phase = 2.0_dp*pi*real(global_y,dp)/real(ny,dp)
      do ii = 1, nx_loc
        ux_i(ii,jj) = amplitude*sin(phase)
      end do
    end do
  end subroutine initialize_shear_wave


  subroutine initialize_density_wave(rho_i, ux_i, uy_i, amplitude)
    real(dp), intent(out) :: rho_i(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(out) :: ux_i(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(out) :: uy_i(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(in) :: amplitude
    integer :: ii, jj, global_x
    real(dp) :: phase_x
    real(dp), parameter :: pi = acos(-1.0_dp)

    ! As with the shear wave, construct the perturbation in global rather than
    ! local coordinates so every decomposition describes the same problem.
    call initialize_uniform(rho_i,ux_i,uy_i)
    do jj = 1, ny_loc
      do ii = 1, nx_loc
        global_x = x_offset + ii - 1
        phase_x = 2.0_dp*pi*real(global_x,dp)/real(nx,dp)
        rho_i(ii,jj) = 1.0_dp + amplitude*sin(phase_x)
      end do
    end do
  end subroutine initialize_density_wave


  subroutine initialize_equilibrium(f_i, rho_i, ux_i, uy_i)
    real(dp), intent(out) :: f_i(q,0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(in) :: rho_i(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(in) :: ux_i(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(in) :: uy_i(0:nx_alloc+1,0:ny_alloc+1)
    integer :: ii, jj, qq
    real(dp) :: cu, u2

    ! Populate each owned node with the second-order, low-Mach D2Q9 equilibrium.
    ! Ghost layers are initially zero and are filled from neighbours before the
    ! first streaming operation.
    f_i = 0.0_dp
    do jj = 1, ny_loc
      do ii = 1, nx_loc
        u2 = ux_i(ii,jj)**2 + uy_i(ii,jj)**2
        do qq = 1, q
          cu = real(cx(qq),dp)*ux_i(ii,jj) + real(cy(qq),dp)*uy_i(ii,jj)
          f_i(qq,ii,jj) = w(qq)*rho_i(ii,jj) * &
            (1.0_dp + cu/cs2 + cu**2/(2.0_dp*cs2**2) - u2/(2.0_dp*cs2))
        end do
      end do
    end do
  end subroutine initialize_equilibrium


  subroutine macro_val(f_m, rho_m, ux_m, uy_m, gx, gy)
    real(dp), intent(in) :: f_m(q,0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(out) :: rho_m(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(out) :: ux_m(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(out) :: uy_m(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(in) :: gx, gy
    integer :: ii, jj, qq
    real(dp) :: momentum_x, momentum_y

    ! Recover density and velocity moments from the particle populations:
    !   rho = sum_q f_q,  rho*u = sum_q f_q*c_q.
    ! The half-force velocity correction is paired with the Guo forcing term in
    ! collide, giving a time-centred body-force discretization.
    rho_m = 0.0_dp
    ux_m = 0.0_dp
    uy_m = 0.0_dp
    do jj = 1, ny_loc
      do ii = 1, nx_loc
        rho_m(ii,jj) = sum(f_m(:,ii,jj))
        momentum_x = 0.0_dp
        momentum_y = 0.0_dp
        do qq = 1, q
          momentum_x = momentum_x + f_m(qq,ii,jj)*real(cx(qq),dp)
          momentum_y = momentum_y + f_m(qq,ii,jj)*real(cy(qq),dp)
        end do
        if (rho_m(ii,jj) > 0.0_dp) then
          ux_m(ii,jj) = momentum_x/rho_m(ii,jj) + 0.5_dp*gx
          uy_m(ii,jj) = momentum_y/rho_m(ii,jj) + 0.5_dp*gy
        end if
      end do
    end do
  end subroutine macro_val


  subroutine collide(f_c, rho_c, ux_c, uy_c, gx, gy)
    real(dp), intent(inout) :: f_c(q,0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(in) :: rho_c(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(in) :: ux_c(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(in) :: uy_c(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(in) :: gx, gy
    integer :: ii, jj, qq
    real(dp) :: feq, source, cu, u2, ci_dot_force
    real(dp) :: force_density_x, force_density_y

    ! Collision is entirely local: relax each population toward equilibrium and
    ! add the Guo source term for acceleration (gx,gy).  Communication occurs
    ! only after all owned cells have been collided.
    do jj = 1, ny_loc
      do ii = 1, nx_loc
        u2 = ux_c(ii,jj)**2 + uy_c(ii,jj)**2
        force_density_x = rho_c(ii,jj)*gx
        force_density_y = rho_c(ii,jj)*gy
        do qq = 1, q
          cu = real(cx(qq),dp)*ux_c(ii,jj) + real(cy(qq),dp)*uy_c(ii,jj)
          feq = w(qq)*rho_c(ii,jj) * &
            (1.0_dp + cu/cs2 + cu**2/(2.0_dp*cs2**2) - u2/(2.0_dp*cs2))
          ci_dot_force = real(cx(qq),dp)*force_density_x + &
            real(cy(qq),dp)*force_density_y
          source = w(qq) * ( &
            (real(cx(qq),dp)-ux_c(ii,jj))*force_density_x/cs2 + &
            (real(cy(qq),dp)-uy_c(ii,jj))*force_density_y/cs2 + &
            cu*ci_dot_force/(cs2**2))
          f_c(qq,ii,jj) = f_c(qq,ii,jj) - &
            relaxation_omega*(f_c(qq,ii,jj)-feq) + &
            (1.0_dp-0.5_dp*relaxation_omega)*source
        end do
      end do
    end do
  end subroutine collide


  subroutine stream(f_s, f_next)
    real(dp), intent(in) :: f_s(q,0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(out) :: f_next(q,0:nx_alloc+1,0:ny_alloc+1)
    integer :: ii, jj, qq

    ! Pull streaming asks where each population at (i,j) came from.  For cells
    ! next to an image boundary, the source index is in the ghost layer already
    ! filled by communicate_ghost_cells.  At a physical wall, those incoming
    ! slots remain zero temporarily and apply_boundaries replaces them below.
    !
    ! Clearing all of f_next also keeps ghost/padding storage deterministic; only
    ! the owned physical region is populated by this loop.
    f_next = 0.0_dp
    do jj = 1, ny_loc
      do ii = 1, nx_loc
        do qq = 1, q
          f_next(qq,ii,jj) = f_s(qq,ii-cx(qq),jj-cy(qq))
        end do
      end do
    end do
  end subroutine stream


  subroutine apply_boundaries(f_post, f_next, u_wall)
    real(dp), intent(in) :: f_post(q,0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(inout) :: f_next(q,0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(in) :: u_wall

    ! Boundary conditions are applied only by images touching a physical global
    ! edge.  Internal interfaces and periodic global edges have nonzero neighbour
    ! ids and were already handled by ghost exchange plus ordinary streaming.
    !
    ! f_post is the post-collision, pre-streaming state.  f_next is the streamed
    ! state.  At a wall, pull streaming could not obtain populations entering the
    ! domain, so the routines below reconstruct exactly those missing directions
    ! from their outgoing opposites in f_post (half-way bounce-back).
    if (left_img == 0 .and. bc_left == BC_WALL) then
      call apply_left_wall(f_post,f_next)
    end if
    if (right_img == 0 .and. bc_right == BC_WALL) then
      call apply_right_wall(f_post,f_next)
    end if
    if (down_img == 0 .and. bc_bottom == BC_WALL) then
      call apply_bottom_wall(f_post,f_next)
    end if
    if (up_img == 0) then
      if (bc_top == BC_WALL) call apply_top_wall(f_post,f_next)
      if (bc_top == BC_MOVING_WALL) then
        call apply_moving_top_wall(f_post,f_next,u_wall)
      end if
    end if
  end subroutine apply_boundaries


  subroutine apply_left_wall(f_post,f_next)
    real(dp), intent(in) :: f_post(q,0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(inout) :: f_next(q,0:nx_alloc+1,0:ny_alloc+1)
    integer :: jj
    ! At x=1, populations 2(E), 6(NE), and 9(SE) would enter from outside.
    ! Reflect the outgoing W, SW, and NW populations respectively.  The wall is
    ! located halfway between the boundary-node centre and the exterior node.
    do jj = 1, ny_loc
      f_next(2,1,jj) = f_post(4,1,jj)
      f_next(6,1,jj) = f_post(8,1,jj)
      f_next(9,1,jj) = f_post(7,1,jj)
    end do
  end subroutine apply_left_wall


  subroutine apply_right_wall(f_post,f_next)
    real(dp), intent(in) :: f_post(q,0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(inout) :: f_next(q,0:nx_alloc+1,0:ny_alloc+1)
    integer :: jj
    ! At x=nx_loc, reconstruct incoming W, NW, and SW populations by reflecting
    ! the corresponding outgoing E, SE, and NE populations.
    do jj = 1, ny_loc
      f_next(4,nx_loc,jj) = f_post(2,nx_loc,jj)
      f_next(7,nx_loc,jj) = f_post(9,nx_loc,jj)
      f_next(8,nx_loc,jj) = f_post(6,nx_loc,jj)
    end do
  end subroutine apply_right_wall


  subroutine apply_bottom_wall(f_post,f_next)
    real(dp), intent(in) :: f_post(q,0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(inout) :: f_next(q,0:nx_alloc+1,0:ny_alloc+1)
    integer :: ii
    ! At y=1, reconstruct the incoming N, NE, and NW populations from the
    ! outgoing S, SW, and SE populations.
    do ii = 1, nx_loc
      f_next(3,ii,1) = f_post(5,ii,1)
      f_next(6,ii,1) = f_post(8,ii,1)
      f_next(7,ii,1) = f_post(9,ii,1)
    end do
  end subroutine apply_bottom_wall


  subroutine apply_top_wall(f_post,f_next)
    real(dp), intent(in) :: f_post(q,0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(inout) :: f_next(q,0:nx_alloc+1,0:ny_alloc+1)
    integer :: ii
    ! At y=ny_loc, reconstruct incoming S, SW, and SE from outgoing N, NE, and
    ! NW populations.  This is no-slip half-way bounce-back for a fixed wall.
    do ii = 1, nx_loc
      f_next(5,ii,ny_loc) = f_post(3,ii,ny_loc)
      f_next(8,ii,ny_loc) = f_post(6,ii,ny_loc)
      f_next(9,ii,ny_loc) = f_post(7,ii,ny_loc)
    end do
  end subroutine apply_top_wall


  subroutine apply_moving_top_wall(f_post,f_next,u_wall)
    real(dp), intent(in) :: f_post(q,0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(inout) :: f_next(q,0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(in) :: u_wall
    integer :: ii
    real(dp) :: rho_wall

    ! Moving-wall bounce-back adds the tangential momentum of a wall travelling
    ! at (+u_wall,0).  The normal population (S) is simply reflected; the two
    ! diagonals receive equal and opposite corrections.  For D2Q9 with cs2=1/3,
    ! the diagonal correction reduces to rho*u_wall/6.
    !
    ! rho_wall is estimated from the local post-collision populations.  At the
    ! top corners, this routine runs after the left/right wall routines and thus
    ! owns the final values of the two populations entering from the top.  This
    ! is the chosen simple corner convention for the moving-lid case.
    do ii = 1, nx_loc
      rho_wall = sum(f_post(:,ii,ny_loc))
      f_next(5,ii,ny_loc) = f_post(3,ii,ny_loc)
      f_next(8,ii,ny_loc) = f_post(6,ii,ny_loc) - rho_wall*u_wall/6.0_dp
      f_next(9,ii,ny_loc) = f_post(7,ii,ny_loc) + rho_wall*u_wall/6.0_dp
    end do
  end subroutine apply_moving_top_wall


  subroutine velocity_residual(ux_c,uy_c,ux_old,uy_old,value)
    real(dp), intent(in) :: ux_c(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(in) :: uy_c(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(inout) :: ux_old(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(inout) :: uy_old(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(out) :: value
    real(dp) :: numerator, denominator

    ! Form a decomposition-independent global relative L2 change in velocity.
    ! Each image contributes only owned cells; co_sum combines those partial
    ! sums on every image before the common convergence decision is made.
    numerator = sum((ux_c(1:nx_loc,1:ny_loc)- &
      ux_old(1:nx_loc,1:ny_loc))**2) + &
      sum((uy_c(1:nx_loc,1:ny_loc)-uy_old(1:nx_loc,1:ny_loc))**2)
    denominator = sum(ux_c(1:nx_loc,1:ny_loc)**2) + &
      sum(uy_c(1:nx_loc,1:ny_loc)**2)
    call co_sum(numerator)
    call co_sum(denominator)
    value = sqrt(numerator/max(denominator,tiny(1.0_dp)))
    ux_old = ux_c
    uy_old = uy_c
  end subroutine velocity_residual

end module lbm_solver
