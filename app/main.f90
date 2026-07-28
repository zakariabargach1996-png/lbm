program main
  use lbm_params
  use domain_decomposition
  use lbm_solver
  use lbm_output
  implicit none

  ! User-facing validation/flow cases.  The numeric values are also passed to
  ! output routines that select the corresponding analytical reference profile.
  integer, parameter :: CASE_SHEAR = 1
  integer, parameter :: CASE_DENSITY = 2
  integer, parameter :: CASE_COUETTE = 3
  integer, parameter :: CASE_POISEUILLE = 4
  integer, parameter :: CASE_LID = 5

  real(dp), parameter :: shear_u0 = 0.05_dp
  real(dp), parameter :: density_amplitude = 0.01_dp
  real(dp), parameter :: target_poiseuille_umax = 0.05_dp
  real(dp), parameter :: convergence_tolerance = 1.0e-8_dp
  integer, parameter :: check_every = 500
  integer, parameter :: wave_output_every = 100
  integer, parameter :: density_output_every = 10

  ! f is the only coarray field: each image exposes its populations so neighbours
  ! can read edge cells into their ghost layers.  f_new and all macroscopic
  ! fields are private to an image.  Every array is indexed with one ghost layer
  ! on each side; owned cells are always 1:nx_loc by 1:ny_loc.
  real(dp), allocatable :: f(:,:,:)[:]
  real(dp), allocatable :: f_new(:,:,:)
  real(dp), allocatable :: rho(:,:), ux(:,:), uy(:,:)
  real(dp), allocatable :: ux_previous(:,:), uy_previous(:,:)

  integer :: case_id, nsteps, minimum_steps, step, final_step
  integer :: profile_output_interval, field_output_interval
  integer :: clock_start, clock_stop, clock_rate
  real(dp) :: force_x, force_y, residual, wall_seconds, mlups
  real(dp) :: relaxation_omega, viscosity, wall_velocity
  logical :: converged, fixed_step_run

  ! Boundary types must be configured before setup_images, because the latter
  ! uses them to decide whether global-edge neighbours wrap or become walls.
  call choose_case(case_id)
  call configure_runtime(relaxation_omega, wall_velocity)
  viscosity = cs2*(1.0_dp/relaxation_omega - 0.5_dp)
  call configure_solver(relaxation_omega)
  call configure_output(relaxation_omega, viscosity)
  call configure_case(case_id, force_x, force_y)
  call setup_images()
  call calculate_run_control(case_id, nsteps, minimum_steps)
  call apply_step_controls(nsteps, minimum_steps, fixed_step_run)
  profile_output_interval = max(check_every,max(1,minimum_steps/4))
  ! The steady solver often converges well before its conservative maximum.
  ! Sampling every half diffusion time gives several useful lid snapshots.
  field_output_interval = max(check_every,max(1,minimum_steps/2))

  ! Coarray bounds must agree on all images.  nx_alloc/ny_alloc are therefore
  ! the maximum block dimensions; nx_loc/ny_loc delimit the part actually owned
  ! by this image when an uneven decomposition leaves smaller edge blocks.
  allocate(f(q,0:nx_alloc+1,0:ny_alloc+1)[*])
  allocate(f_new(q,0:nx_alloc+1,0:ny_alloc+1))
  allocate(rho(0:nx_alloc+1,0:ny_alloc+1))
  allocate(ux(0:nx_alloc+1,0:ny_alloc+1))
  allocate(uy(0:nx_alloc+1,0:ny_alloc+1))
  allocate(ux_previous(0:nx_alloc+1,0:ny_alloc+1))
  allocate(uy_previous(0:nx_alloc+1,0:ny_alloc+1))

#ifndef LBM_PERFORMANCE_ONLY
  select case (case_id)
  case (CASE_SHEAR)
    call initialize_shear_wave(rho, ux, uy, shear_u0)
  case (CASE_DENSITY)
    call initialize_density_wave(rho, ux, uy, density_amplitude)
  case default
    call initialize_uniform(rho, ux, uy)
  end select

  call initialize_equilibrium(f, rho, ux, uy)
  f_new = 0.0_dp
  ux_previous = ux
  uy_previous = uy

  ! Seed halos before the first time step.  This is needed because the first
  ! streaming operation may pull populations across an image/periodic boundary.
  call communicate_ghost_cells(f)

  if (my_img == 1) then
    print '(a,f8.5,a,es12.4)', "omega: ", relaxation_omega, &
      ", analytical viscosity: ", viscosity
    print '(a,i0)', "Maximum steps: ", nsteps
    if (case_id >= CASE_COUETTE) then
      print '(a,es12.4)', "Convergence tolerance: ", convergence_tolerance
    end if
  end if

  select case (case_id)
  case (CASE_SHEAR)
    call write_shear_wave_decay(f, 0, shear_u0, "shear_wave_decay.txt")
    call write_shear_profile(rho, ux, uy, 0, "shear_profiles.txt")
  case (CASE_DENSITY)
    call write_density_profile(rho, ux, uy, 0, "density_profiles.txt")
  case (CASE_COUETTE)
    call write_channel_history(ux, 0, "couette_evolution.txt", &
      case_id, force_x, wall_velocity)
  case (CASE_POISEUILLE)
    call write_channel_history(ux, 0, "poiseuille_evolution.txt", &
      case_id, force_x, wall_velocity)
  case (CASE_LID)
    call write_flow_snapshot(rho, ux, uy, 0)
  end select
#endif

  converged = .false.
  residual = huge(1.0_dp)
  final_step = 0

  sync all
  call system_clock(clock_start,clock_rate)

  ! One lattice-Boltzmann time step:
  !   1. recover rho and u from the current populations;
  !   2. collide owned cells locally (and apply any body force);
  !   3. exchange post-collision edge data into one-cell ghost layers;
  !   4. pull-stream from local or ghost source cells;
  !   5. reconstruct populations that enter through physical walls;
  !   6. promote the completed buffer to the current state.
  !
  ! Exchanging after collision is essential: streaming must transport the
  ! neighbour's new post-collision populations, not its previous-step values.
  do step = 1, nsteps
    call macro_val(f, rho, ux, uy, force_x, force_y)
    call collide(f, rho, ux, uy, force_x, force_y)
    call communicate_ghost_cells(f)
    call stream(f, f_new)
    call apply_boundaries(f, f_new, wall_velocity)
    f = f_new
    final_step = step

#ifndef LBM_PERFORMANCE_ONLY
    if (.not. fixed_step_run .and. case_id == CASE_SHEAR .and. &
        mod(step,wave_output_every) == 0) then
      call macro_val(f, rho, ux, uy, force_x, force_y)
      call write_shear_wave_decay(f, step, shear_u0, "shear_wave_decay.txt")
      call write_shear_profile(rho, ux, uy, step, "shear_profiles.txt")
    end if

    if (.not. fixed_step_run .and. case_id == CASE_DENSITY .and. &
        mod(step,density_output_every) == 0) then
      call macro_val(f, rho, ux, uy, force_x, force_y)
      call write_density_profile(rho, ux, uy, step, "density_profiles.txt")
    end if

    if (.not. fixed_step_run .and. &
        (case_id == CASE_COUETTE .or. case_id == CASE_POISEUILLE) .and. &
        mod(step,profile_output_interval) == 0) then
      call macro_val(f, rho, ux, uy, force_x, force_y)
      if (case_id == CASE_COUETTE) then
        call write_channel_history(ux, step, "couette_evolution.txt", &
          case_id, force_x, wall_velocity)
      else
        call write_channel_history(ux, step, "poiseuille_evolution.txt", &
          case_id, force_x, wall_velocity)
      end if
    end if

    if (.not. fixed_step_run .and. case_id == CASE_LID .and. &
        mod(step,field_output_interval) == 0) then
      call macro_val(f, rho, ux, uy, force_x, force_y)
      call write_flow_snapshot(rho, ux, uy, step)
    end if
#endif

    if (.not. fixed_step_run .and. case_id >= CASE_COUETTE .and. &
        mod(step,check_every) == 0) then
      call macro_val(f, rho, ux, uy, force_x, force_y)
      call velocity_residual(ux, uy, ux_previous, uy_previous, residual)
      if (my_img == 1 .and. mod(step,10*check_every) == 0) then
        print '(a,i0,a,es12.4)', "step ", step, ", residual ", residual
      end if
      if (step >= minimum_steps .and. residual < convergence_tolerance) then
        converged = .true.
        exit
      end if
    end if
  end do

  sync all
  call system_clock(clock_stop)
  ! The slowest image determines parallel wall time.  co_max sends that value to
  ! image 1, where global throughput is reported in million lattice updates/s.
  wall_seconds = real(clock_stop-clock_start,dp)/real(clock_rate,dp)
  call co_max(wall_seconds,result_image=1)
  if (my_img == 1) then
    mlups = real(nx,dp)*real(ny,dp)*real(final_step,dp) / &
      (1.0e6_dp*max(wall_seconds,tiny(1.0_dp)))
    call report_performance(case_id,final_step,wall_seconds,mlups)
  end if

#ifndef LBM_PERFORMANCE_ONLY
  call macro_val(f, rho, ux, uy, force_x, force_y)

  select case (case_id)
  case (CASE_SHEAR)
    if (fixed_step_run .or. mod(final_step,wave_output_every) /= 0) then
      call write_shear_wave_decay(f, final_step, shear_u0, "shear_wave_decay.txt")
      call write_shear_profile(rho, ux, uy, final_step, "shear_profiles.txt")
    end if
    call write_viscosity_measurement(f, final_step, shear_u0)
  case (CASE_DENSITY)
    if (fixed_step_run .or. mod(final_step,density_output_every) /= 0) then
      call write_density_profile(rho, ux, uy, final_step, "density_profiles.txt")
    end if
  case (CASE_COUETTE)
    if (fixed_step_run .or. mod(final_step,profile_output_interval) /= 0) then
      call write_channel_history(ux, final_step, "couette_evolution.txt", &
        case_id, force_x, wall_velocity)
    end if
  case (CASE_POISEUILLE)
    if (fixed_step_run .or. mod(final_step,profile_output_interval) /= 0) then
      call write_channel_history(ux, final_step, "poiseuille_evolution.txt", &
        case_id, force_x, wall_velocity)
    end if
  case (CASE_LID)
    if (fixed_step_run .or. mod(final_step,field_output_interval) /= 0) then
      call write_flow_snapshot(rho, ux, uy, final_step)
    end if
  end select

  select case (case_id)
  case (CASE_COUETTE)
    call write_couette_profile(ux, wall_velocity, "couette_profile.txt")
  case (CASE_POISEUILLE)
    call write_poiseuille_profile(ux, force_x, "poiseuille_profile.txt")
  case (CASE_LID)
    call write_flow_field(rho, ux, uy, "moving_lid_field.txt")
  end select
#endif

  if (my_img == 1) then
    print '(a,i0)', "Completed steps: ", final_step
    if (case_id >= CASE_COUETTE) then
      if (.not. fixed_step_run .and. final_step >= check_every) then
        print '(a,es12.4)', "Final residual: ", residual
        if (converged) then
          print '(a)', "Convergence criterion reached."
        else
          print '(a)', "Maximum step limit reached before convergence."
        end if
      else
        print '(a)', "Residual not sampled during this short fixed-step run."
      end if
    end if
  end if

contains

  subroutine choose_case(selected)
    integer, intent(out) :: selected
    character(len=64) :: argument
    integer :: status

    ! Prefer a command-line case name/number.  In interactive use, only image 1
    ! reads stdin; co_broadcast gives every image the identical selection.
    selected = 0
    argument = ""
    call get_command_argument(1, argument, status=status)

    if (status == 0 .and. len_trim(argument) > 0) then
      call lowercase(argument)
      select case (trim(argument))
      case ("1", "shear", "shear-wave", "shear_wave")
        selected = CASE_SHEAR
      case ("2", "density", "density-wave", "density_wave")
        selected = CASE_DENSITY
      case ("3", "couette")
        selected = CASE_COUETTE
      case ("4", "poiseuille")
        selected = CASE_POISEUILLE
      case ("5", "lid", "moving-lid", "moving_lid")
        selected = CASE_LID
      end select
    end if

    if (selected == 0 .and. this_image() == 1) then
      print '(a)', "Select an LBM case:"
      print '(a)', "  1 - Shear-wave decay"
      print '(a)', "  2 - Density-wave animation"
      print '(a)', "  3 - Couette flow"
      print '(a)', "  4 - Poiseuille flow"
      print '(a)', "  5 - Moving-lid cavity"
      read (*,*,iostat=status) selected
      if (status /= 0) selected = 0
    end if

    call co_broadcast(selected,source_image=1)
    if (selected < CASE_SHEAR .or. selected > CASE_LID) then
      if (this_image() == 1) print '(a)', "Invalid case selection."
      error stop
    end if
  end subroutine choose_case


  subroutine configure_runtime(omega_value, u_wall)
    real(dp), intent(out) :: omega_value, u_wall
    character(len=64) :: argument
    integer :: io_status

    ! Parse optional parameters on image 1 only, then broadcast them.  Keeping
    ! command-line I/O on one image avoids inconsistent runtime configuration.
    omega_value = default_omega
    u_wall = 0.05_dp

    if (this_image() == 1) then
      call get_command_argument(2, argument)
      if (len_trim(argument) > 0) then
        read(argument,*,iostat=io_status) omega_value
        if (io_status /= 0) omega_value = -1.0_dp
      end if

      call get_command_argument(3, argument)
      if (len_trim(argument) > 0) then
        read(argument,*,iostat=io_status) u_wall
        if (io_status /= 0) omega_value = -1.0_dp
      end if
    end if

    call co_broadcast(omega_value,source_image=1)
    call co_broadcast(u_wall,source_image=1)
    if (omega_value <= 0.0_dp .or. omega_value >= 2.0_dp) then
      if (this_image() == 1) then
        print '(a)', "omega must satisfy 0 < omega < 2."
        print '(a)', &
          "Usage: milestone4 CASE [OMEGA] [WALL_VELOCITY] [MAX_STEPS]"
      end if
      error stop
    end if
  end subroutine configure_runtime


  subroutine lowercase(text)
    character(len=*), intent(inout) :: text
    integer :: ii, code

    do ii = 1, len_trim(text)
      code = iachar(text(ii:ii))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        text(ii:ii) = achar(code + iachar('a') - iachar('A'))
      end if
    end do
  end subroutine lowercase


  subroutine apply_step_controls(maximum_steps,min_steps,is_fixed)
    integer, intent(inout) :: maximum_steps, min_steps
    logical, intent(out) :: is_fixed
    character(len=64) :: argument
    integer :: requested_steps, step_limit, io_status, env_status

    ! Supplying MAX_STEPS selects deterministic benchmark mode.  Diagnostics and
    ! early convergence are disabled so timings cover exactly the requested work.
    requested_steps = 0
    is_fixed = .false.
    if (this_image() == 1) then
      call get_command_argument(4,argument)
      if (len_trim(argument) > 0) then
        read(argument,*,iostat=io_status) requested_steps
        if (io_status /= 0 .or. requested_steps < 1) requested_steps = -1
      end if
    end if
    call co_broadcast(requested_steps,source_image=1)
    if (requested_steps < 0) then
      if (this_image() == 1) print '(a)', "MAX_STEPS must be a positive integer."
      error stop
    end if
    if (requested_steps > 0) then
      maximum_steps = requested_steps
      is_fixed = .true.
      ! Fixed-step benchmarks suppress in-loop diagnostics and convergence.
      min_steps = requested_steps + 1
      if (this_image() == 1) print '(a)', &
        "Fixed-step timing mode: in-loop output and convergence checks disabled."
    end if

    ! LBM_STEP_LIMIT is a safety cap for physical runs. Unlike the positional
    ! MAX_STEPS benchmark argument, it preserves diagnostics and early stopping.
    step_limit = 0
    if (this_image() == 1 .and. .not. is_fixed) then
      argument = ""
      call get_environment_variable("LBM_STEP_LIMIT",argument,status=env_status)
      if (env_status == 0 .and. len_trim(argument) > 0) then
        read(argument,*,iostat=io_status) step_limit
        if (io_status /= 0 .or. step_limit < 1) step_limit = -1
      end if
    end if
    call co_broadcast(step_limit,source_image=1)
    if (step_limit < 0) then
      if (this_image() == 1) then
        print '(a)', "LBM_STEP_LIMIT must be a positive integer."
      end if
      error stop
    end if
    if (step_limit > 0) then
      maximum_steps = min(maximum_steps,step_limit)
      if (this_image() == 1) then
        print '(a,i0)', "Convergence-aware step limit: ", maximum_steps
      end if
    end if
  end subroutine apply_step_controls


  subroutine report_performance(selected,completed_steps,seconds,rate_mlups)
    integer, intent(in) :: selected, completed_steps
    real(dp), intent(in) :: seconds, rate_mlups
    integer :: io_unit
    logical :: exists
    character(len=16) :: case_name

    select case (selected)
    case (CASE_SHEAR);      case_name = "shear"
    case (CASE_DENSITY);    case_name = "density"
    case (CASE_COUETTE);    case_name = "couette"
    case (CASE_POISEUILLE); case_name = "poiseuille"
    case (CASE_LID);        case_name = "lid"
    end select

    print '(a,f12.6)', "Wall-clock time [s]: ", seconds
    print '(a,f12.3)', "Performance [MLUPS]: ", rate_mlups
    inquire(file="performance.csv",exist=exists)
    if (exists) then
      open(newunit=io_unit,file="performance.csv",status="old", &
        position="append",action="write")
    else
      open(newunit=io_unit,file="performance.csv",status="new",action="write")
      write(io_unit,'(a)') "case,nx,ny,images,px,py,steps,wall_seconds,mlups"
    end if
    write(io_unit,'(a,6(",",i0),2(",",es24.16))') trim(case_name),nx,ny, &
      nimgs,px,py,completed_steps,seconds,rate_mlups
    close(io_unit)
  end subroutine report_performance


  subroutine configure_case(selected, gx, gy)
    integer, intent(in) :: selected
    real(dp), intent(out) :: gx, gy

    ! These are physical global-edge conditions, not per-image conditions.
    ! setup_images later translates them into neighbour ids for each block:
    ! periodic edges wrap during halo exchange; walls receive neighbour id zero
    ! and are reconstructed after streaming.
    gx = 0.0_dp
    gy = 0.0_dp

    select case (selected)
    case (CASE_SHEAR, CASE_DENSITY)
      bc_left = BC_PERIODIC
      bc_right = BC_PERIODIC
      bc_bottom = BC_PERIODIC
      bc_top = BC_PERIODIC
    case (CASE_COUETTE)
      bc_left = BC_PERIODIC
      bc_right = BC_PERIODIC
      bc_bottom = BC_WALL
      bc_top = BC_MOVING_WALL
    case (CASE_POISEUILLE)
      bc_left = BC_PERIODIC
      bc_right = BC_PERIODIC
      bc_bottom = BC_WALL
      bc_top = BC_WALL
      ! Plane Poiseuille flow has umax = gx*H^2/(8*nu).
      gx = 8.0_dp*viscosity*target_poiseuille_umax/real(ny,dp)**2
    case (CASE_LID)
      bc_left = BC_WALL
      bc_right = BC_WALL
      bc_bottom = BC_WALL
      bc_top = BC_MOVING_WALL
    end select

    if (this_image() == 1) then
      select case (selected)
      case (CASE_SHEAR)
        print '(a)', "Case: shear-wave decay"
      case (CASE_DENSITY)
        print '(a)', "Case: density-wave animation"
      case (CASE_COUETTE)
        print '(a)', "Case: Couette flow"
      case (CASE_POISEUILLE)
        print '(a,es12.4)', "Case: Poiseuille flow, acceleration = ", gx
      case (CASE_LID)
        print '(a,f10.3)', "Case: moving-lid cavity, Reynolds number = ", &
          wall_velocity*real(nx,dp)/viscosity
      end select
    end if
  end subroutine configure_case


  subroutine calculate_run_control(selected, maximum_steps, min_steps)
    integer, intent(in) :: selected
    integer, intent(out) :: maximum_steps, min_steps
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: slow_time, wave_number

    ! Estimate the slowest viscous diffusion time.  Steady cases must run for at
    ! least one such time before the residual is allowed to stop the simulation.
    slow_time = real(max(nx,ny),dp)**2/(pi*pi*viscosity)
    wave_number = 2.0_dp*pi/real(ny,dp)
    min_steps = 0

    select case (selected)
    case (CASE_SHEAR)
      maximum_steps = ceiling(-log(0.2_dp)/(viscosity*wave_number**2))
    case (CASE_DENSITY)
      ! Five acoustic periods provide enough frames for a useful animation.
      maximum_steps = ceiling(5.0_dp*real(nx,dp)/sqrt(cs2))
    case default
      min_steps = ceiling(slow_time)
      maximum_steps = ceiling(-log(convergence_tolerance)*slow_time)
    end select
  end subroutine calculate_run_control





end program main
