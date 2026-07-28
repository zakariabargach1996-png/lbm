module lbm_output
  use lbm_params
  use domain_decomposition
  implicit none
  private

  public :: configure_output
  public :: write_shear_wave_decay, write_viscosity_measurement
  public :: write_shear_profile, write_density_profile, write_channel_history
  public :: write_couette_profile, write_poiseuille_profile
  public :: write_flow_field, write_flow_snapshot

  integer, parameter :: CASE_COUETTE = 3
  ! Cached physical parameters used to construct analytical reference solutions.
  real(dp) :: relaxation_omega, viscosity

contains

  subroutine configure_output(omega,nu)
    real(dp), intent(in) :: omega, nu

    relaxation_omega = omega
    viscosity = nu
  end subroutine configure_output
  subroutine write_shear_wave_decay(f_s,current_step,u0,filename)
    real(dp), intent(in) :: f_s(q,0:nx_alloc+1,0:ny_alloc+1)
    integer, intent(in) :: current_step
    real(dp), intent(in) :: u0
    character(len=*), intent(in) :: filename
    integer :: io_unit
    real(dp) :: amplitude, theoretical, relative_error
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: wave_number = 2.0_dp*pi/real(ny,dp)

    call measure_shear_amplitude(f_s,amplitude)

    if (my_img == 1) then
      theoretical = u0*exp(-viscosity*wave_number**2*real(current_step,dp))
      relative_error = abs(amplitude-theoretical)/max(abs(theoretical),tiny(1.0_dp))
      if (current_step == 0) then
        open(newunit=io_unit,file=filename,status="replace",action="write")
        write(io_unit,'(a)') "# step amplitude theoretical relative_error"
      else
        open(newunit=io_unit,file=filename,status="old",position="append",action="write")
      end if
      write(io_unit,'(i10,3(1x,es24.16))') current_step, amplitude, &
        theoretical, relative_error
      close(io_unit)
    end if
  end subroutine write_shear_wave_decay


  subroutine measure_shear_amplitude(f_s,amplitude)
    real(dp), intent(in) :: f_s(q,0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(out) :: amplitude
    integer :: ii, jj, qq, global_y
    real(dp) :: rho_cell, velocity, phase, projection
    real(dp), parameter :: pi = acos(-1.0_dp)

    ! Project ux onto the initialized sine mode.  Every image accumulates its
    ! owned cells using global coordinates; the reduction assembles the global
    ! amplitude on image 1 without first gathering the complete field.
    projection = 0.0_dp
    do jj = 1, ny_loc
      global_y = y_offset + jj - 1
      phase = 2.0_dp*pi*real(global_y,dp)/real(ny,dp)
      do ii = 1, nx_loc
        rho_cell = sum(f_s(:,ii,jj))
        velocity = 0.0_dp
        do qq = 1, q
          velocity = velocity + f_s(qq,ii,jj)*real(cx(qq),dp)
        end do
        projection = projection + velocity*sin(phase)/rho_cell
      end do
    end do
    call co_sum(projection,result_image=1)
    amplitude = 0.0_dp
    if (my_img == 1) amplitude = 2.0_dp*projection/real(nx*ny,dp)
  end subroutine measure_shear_amplitude


  subroutine write_viscosity_measurement(f_s,current_step,u0)
    real(dp), intent(in) :: f_s(q,0:nx_alloc+1,0:ny_alloc+1)
    integer, intent(in) :: current_step
    real(dp), intent(in) :: u0
    real(dp) :: amplitude, measured_viscosity, wave_number
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: io_unit
    character(len=64) :: filename

    call measure_shear_amplitude(f_s,amplitude)
    if (my_img == 1) then
      wave_number = 2.0_dp*pi/real(ny,dp)
      measured_viscosity = -log(abs(amplitude/u0))/ &
        (wave_number**2*real(current_step,dp))
      write(filename,'("viscosity_omega_",f5.3,".txt")') relaxation_omega
      open(newunit=io_unit,file=trim(filename),status="replace",action="write")
      write(io_unit,'(a)') "# omega measured_viscosity analytical_viscosity relative_error"
      write(io_unit,'(4(es24.16,1x))') relaxation_omega, measured_viscosity, &
        viscosity, abs(measured_viscosity-viscosity)/viscosity
      close(io_unit)
      print '(a,es12.4)', "Measured viscosity: ", measured_viscosity
    end if
  end subroutine write_viscosity_measurement


  subroutine write_shear_profile(rho_c,ux_c,uy_c,current_step,filename)
    real(dp), intent(in) :: rho_c(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(in) :: ux_c(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(in) :: uy_c(0:nx_alloc+1,0:ny_alloc+1)
    integer, intent(in) :: current_step
    character(len=*), intent(in) :: filename
    real(dp) :: rho_profile(ny), ux_profile(ny), uy_profile(ny)
    integer :: jj, global_y, io_unit

    ! Each global y entry receives contributions only from images owning that
    ! row.  co_sum combines the x-block partial sums into a global mean profile.
    rho_profile = 0.0_dp
    ux_profile = 0.0_dp
    uy_profile = 0.0_dp
    do jj = 1, ny_loc
      global_y = y_offset + jj
      rho_profile(global_y) = sum(rho_c(1:nx_loc,jj))
      ux_profile(global_y) = sum(ux_c(1:nx_loc,jj))
      uy_profile(global_y) = sum(uy_c(1:nx_loc,jj))
    end do
    call co_sum(rho_profile,result_image=1)
    call co_sum(ux_profile,result_image=1)
    call co_sum(uy_profile,result_image=1)

    if (my_img == 1) then
      if (current_step == 0) then
        open(newunit=io_unit,file=filename,status="replace",action="write")
        write(io_unit,'(a)') "# step y mean_rho mean_ux mean_uy"
      else
        open(newunit=io_unit,file=filename,status="old",position="append",action="write")
      end if
      do jj = 1, ny
        write(io_unit,'(i10,1x,f12.4,3(1x,es24.16))') current_step, &
          real(jj-1,dp),rho_profile(jj)/real(nx,dp), &
          ux_profile(jj)/real(nx,dp),uy_profile(jj)/real(nx,dp)
      end do
      write(io_unit,*)
      close(io_unit)
    end if
  end subroutine write_shear_profile


  subroutine write_density_profile(rho_c,ux_c,uy_c,current_step,filename)
    real(dp), intent(in) :: rho_c(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(in) :: ux_c(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(in) :: uy_c(0:nx_alloc+1,0:ny_alloc+1)
    integer, intent(in) :: current_step
    character(len=*), intent(in) :: filename
    real(dp) :: rho_profile(nx), ux_profile(nx), uy_profile(nx)
    integer :: ii, global_x, io_unit

    ! Assemble an x profile analogously, reducing partial sums over y blocks.
    rho_profile = 0.0_dp
    ux_profile = 0.0_dp
    uy_profile = 0.0_dp
    do ii = 1, nx_loc
      global_x = x_offset + ii
      rho_profile(global_x) = sum(rho_c(ii,1:ny_loc))
      ux_profile(global_x) = sum(ux_c(ii,1:ny_loc))
      uy_profile(global_x) = sum(uy_c(ii,1:ny_loc))
    end do
    call co_sum(rho_profile,result_image=1)
    call co_sum(ux_profile,result_image=1)
    call co_sum(uy_profile,result_image=1)

    if (my_img == 1) then
      if (current_step == 0) then
        open(newunit=io_unit,file=filename,status="replace",action="write")
        write(io_unit,'(a)') "# step x mean_rho mean_ux mean_uy"
      else
        open(newunit=io_unit,file=filename,status="old",position="append",action="write")
      end if
      do ii = 1, nx
        write(io_unit,'(i10,1x,f12.4,3(1x,es24.16))') current_step, &
          real(ii-1,dp),rho_profile(ii)/real(ny,dp), &
          ux_profile(ii)/real(ny,dp),uy_profile(ii)/real(ny,dp)
      end do
      write(io_unit,*)
      close(io_unit)
    end if
  end subroutine write_density_profile


  subroutine collect_x_profile(ux_c,profile)
    real(dp), intent(in) :: ux_c(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(out) :: profile(ny)
    integer :: jj, global_y

    ! Compute the streamwise mean at each global y coordinate.  The arrays are
    ! replicated temporaries, and only image 1 receives/normalizes the reduction.
    profile = 0.0_dp
    do jj = 1, ny_loc
      global_y = y_offset + jj
      profile(global_y) = sum(ux_c(1:nx_loc,jj))
    end do
    call co_sum(profile,result_image=1)
    if (my_img == 1) profile = profile/real(nx,dp)
  end subroutine collect_x_profile


  subroutine write_channel_history(ux_c,current_step,filename,selected,gx,u_wall)
    real(dp), intent(in) :: ux_c(0:nx_alloc+1,0:ny_alloc+1)
    integer, intent(in) :: current_step, selected
    character(len=*), intent(in) :: filename
    real(dp), intent(in) :: gx, u_wall
    real(dp) :: profile(ny), y_position, theoretical
    integer :: jj, io_unit

    ! Half-way bounce-back places the physical walls half a lattice spacing from
    ! the first/last fluid-node centres.  Thus fluid nodes have y=j-1/2 while the
    ! analytical wall locations written around them are y=0 and y=ny.
    call collect_x_profile(ux_c,profile)
    if (my_img == 1) then
      if (current_step == 0) then
        open(newunit=io_unit,file=filename,status="replace",action="write")
        write(io_unit,'(a)') "# step y numerical_ux steady_analytical_ux"
      else
        open(newunit=io_unit,file=filename,status="old",position="append",action="write")
      end if

      write(io_unit,'(i10,3(1x,es24.16))') current_step,0.0_dp,0.0_dp,0.0_dp
      do jj = 1, ny
        y_position = real(jj,dp)-0.5_dp
        if (selected == CASE_COUETTE) then
          theoretical = u_wall*y_position/real(ny,dp)
        else
          theoretical = gx*y_position*(real(ny,dp)-y_position)/ &
            (2.0_dp*viscosity)
        end if
        write(io_unit,'(i10,3(1x,es24.16))') current_step,y_position, &
          profile(jj),theoretical
      end do
      if (selected == CASE_COUETTE) then
        theoretical = u_wall
      else
        theoretical = 0.0_dp
      end if
      write(io_unit,'(i10,3(1x,es24.16))') current_step,real(ny,dp), &
        theoretical,theoretical
      write(io_unit,*)
      close(io_unit)
    end if
  end subroutine write_channel_history


  subroutine write_couette_profile(ux_c,u_wall,filename)
    real(dp), intent(in) :: ux_c(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(in) :: u_wall
    character(len=*), intent(in) :: filename
    real(dp) :: profile(ny), y_position, theoretical, error
    real(dp) :: l2_numerator, l2_denominator, linf_error
    integer :: jj, io_unit

    call collect_x_profile(ux_c,profile)
    if (my_img == 1) then
      l2_numerator = 0.0_dp
      l2_denominator = 0.0_dp
      linf_error = 0.0_dp
      open(newunit=io_unit,file=filename,status="replace",action="write")
      write(io_unit,'(a)') "# y numerical_ux theoretical_ux absolute_error"
      write(io_unit,'(4(es24.16,1x))') 0.0_dp,0.0_dp,0.0_dp,0.0_dp
      do jj = 1, ny
        y_position = real(jj,dp)-0.5_dp
        theoretical = u_wall*y_position/real(ny,dp)
        error = abs(profile(jj)-theoretical)
        l2_numerator = l2_numerator + (profile(jj)-theoretical)**2
        l2_denominator = l2_denominator + theoretical**2
        linf_error = max(linf_error,error)
        write(io_unit,'(4(es24.16,1x))') y_position,profile(jj),theoretical,error
      end do
      write(io_unit,'(4(es24.16,1x))') real(ny,dp),u_wall,u_wall,0.0_dp
      close(io_unit)
      print '(a,es12.4)', "Couette relative L2 error: ", &
        sqrt(l2_numerator/max(l2_denominator,tiny(1.0_dp)))
      print '(a,es12.4)', "Couette Linf error:       ", linf_error
    end if
  end subroutine write_couette_profile


  subroutine write_poiseuille_profile(ux_c,gx,filename)
    real(dp), intent(in) :: ux_c(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(in) :: gx
    character(len=*), intent(in) :: filename
    real(dp) :: profile(ny), y_position, theoretical, error
    real(dp) :: l2_numerator, l2_denominator, linf_error
    integer :: jj, io_unit

    call collect_x_profile(ux_c,profile)
    if (my_img == 1) then
      l2_numerator = 0.0_dp
      l2_denominator = 0.0_dp
      linf_error = 0.0_dp
      open(newunit=io_unit,file=filename,status="replace",action="write")
      write(io_unit,'(a)') "# y numerical_ux theoretical_ux absolute_error"
      write(io_unit,'(4(es24.16,1x))') 0.0_dp,0.0_dp,0.0_dp,0.0_dp
      do jj = 1, ny
        y_position = real(jj,dp)-0.5_dp
        theoretical = gx*y_position*(real(ny,dp)-y_position)/(2.0_dp*viscosity)
        error = abs(profile(jj)-theoretical)
        l2_numerator = l2_numerator + (profile(jj)-theoretical)**2
        l2_denominator = l2_denominator + theoretical**2
        linf_error = max(linf_error,error)
        write(io_unit,'(4(es24.16,1x))') y_position,profile(jj),theoretical,error
      end do
      write(io_unit,'(4(es24.16,1x))') real(ny,dp),0.0_dp,0.0_dp,0.0_dp
      close(io_unit)
      print '(a,es12.4)', "Poiseuille relative L2 error: ", &
        sqrt(l2_numerator/max(l2_denominator,tiny(1.0_dp)))
      print '(a,es12.4)', "Poiseuille Linf error:       ", linf_error
    end if
  end subroutine write_poiseuille_profile


  subroutine write_flow_field(rho_c,ux_c,uy_c,filename)
    real(dp), intent(in) :: rho_c(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(in) :: ux_c(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(in) :: uy_c(0:nx_alloc+1,0:ny_alloc+1)
    character(len=*), intent(in) :: filename
    real(dp), allocatable :: global_rho(:,:), global_ux(:,:), global_uy(:,:)
    integer :: ii, jj, global_x, global_y, io_unit

    ! Assemble a full field through reductions: each image writes its owned block
    ! into otherwise-zero global work arrays, then co_sum overlays the disjoint
    ! blocks on image 1.  Ghosts and allocation padding are never included.
    allocate(global_rho(nx,ny),global_ux(nx,ny),global_uy(nx,ny))
    global_rho = 0.0_dp
    global_ux = 0.0_dp
    global_uy = 0.0_dp
    do jj = 1, ny_loc
      global_y = y_offset + jj
      do ii = 1, nx_loc
        global_x = x_offset + ii
        global_rho(global_x,global_y) = rho_c(ii,jj)
        global_ux(global_x,global_y) = ux_c(ii,jj)
        global_uy(global_x,global_y) = uy_c(ii,jj)
      end do
    end do
    call co_sum(global_rho,result_image=1)
    call co_sum(global_ux,result_image=1)
    call co_sum(global_uy,result_image=1)

    if (my_img == 1) then
      open(newunit=io_unit,file=filename,status="replace",action="write")
      write(io_unit,'(a)') "# x y rho ux uy speed"
      do jj = 1, ny
        do ii = 1, nx
          write(io_unit,'(2(es16.8,1x),4(es24.16,1x))') &
            real(ii,dp)-0.5_dp,real(jj,dp)-0.5_dp,global_rho(ii,jj), &
            global_ux(ii,jj),global_uy(ii,jj), &
            sqrt(global_ux(ii,jj)**2+global_uy(ii,jj)**2)
        end do
        write(io_unit,*)
      end do
      close(io_unit)
    end if
    deallocate(global_rho,global_ux,global_uy)
  end subroutine write_flow_field


  subroutine write_flow_snapshot(rho_c,ux_c,uy_c,current_step)
    real(dp), intent(in) :: rho_c(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(in) :: ux_c(0:nx_alloc+1,0:ny_alloc+1)
    real(dp), intent(in) :: uy_c(0:nx_alloc+1,0:ny_alloc+1)
    integer, intent(in) :: current_step
    character(len=64) :: filename

    write(filename,'("moving_lid_frame_",i10.10,".txt")') current_step
    call write_flow_field(rho_c,ux_c,uy_c,trim(filename))
  end subroutine write_flow_snapshot

end module lbm_output
