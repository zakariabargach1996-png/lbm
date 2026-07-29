module lbm_solver
  use lbm_params
  use domain_decomposition
  implicit none
  private

  public :: configure_solver
  public :: initialize_uniform, initialize_shear_wave, initialize_density_wave
  public :: initialize_equilibrium, macro_val, collide, stream
  public :: apply_boundaries, velocity_residual

  ! BGK collision frequency configured at startup.
  real(dp) :: relaxation_omega

contains

  subroutine configure_solver(omega)
    real(dp), intent(in) :: omega

    ! Store the relaxation parameter used by the collision operator.
    relaxation_omega = omega
  end subroutine configure_solver
  subroutine initialize_uniform(rho_i, ux_i, uy_i)
    real(dp), intent(out) :: rho_i(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(out) :: ux_i(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(out) :: uy_i(0:nx_alloc+1,0:ny_alloc+1)

    ! Initialize a stationary, unit-density fluid.
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

    ! Initialize a globally continuous sinusoidal shear velocity.
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

    ! Initialize a globally continuous sinusoidal density perturbation.
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

    ! Convert macroscopic fields to equilibrium D2Q9 populations.
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

    ! Recover density and force-corrected velocity from the populations.
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

    ! Apply BGK collision and the Guo body-force source term.
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

    ! Pull post-collision populations from local cells or ghost layers.
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

    ! Apply stationary or moving bounce-back at physical domain edges.
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
    ! Reflect populations entering through the left wall.
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
    ! Reflect populations entering through the right wall.
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
    ! Reflect populations entering through the bottom wall.
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
    ! Reflect populations entering through a stationary top wall.
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

    ! Reflect top-wall populations and add the lid's tangential momentum.
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

    ! Compute the global relative L2 velocity change and update the reference.
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
