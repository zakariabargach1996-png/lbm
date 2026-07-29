module domain_decomposition
  use lbm_params
  implicit none

  ! Decompose the global lattice over a two-dimensional coarray image grid.

  integer :: nimgs
  integer :: my_img

  integer :: px, py
  integer :: ix_img, iy_img

  integer :: nx_loc, ny_loc
  integer :: nx_alloc, ny_alloc
  integer :: x_offset, y_offset

  ! A zero neighbour identifies a physical boundary.
  integer :: left_img, right_img
  integer :: down_img, up_img

contains

  subroutine setup_images()
    implicit none

    ! Build local extents, offsets, and neighbours for this image.
    nimgs  = num_images()
    my_img = this_image()

    call choose_image_grid(nimgs, px, py)

    ix_img = mod(my_img - 1, px) + 1
    iy_img = (my_img - 1) / px + 1

    nx_loc = block_size(nx, px, ix_img)
    ny_loc = block_size(ny, py, iy_img)
    nx_alloc = (nx + px - 1) / px
    ny_alloc = (ny + py - 1) / py
    x_offset = block_offset(nx, px, ix_img)
    y_offset = block_offset(ny, py, iy_img)



   if (ix_img == 1) then
      if (bc_left == BC_PERIODIC) then
         left_img = image_id(px, iy_img)
      else
         left_img = 0
      end if
   else
      left_img = image_id(ix_img - 1, iy_img)
   end if


   if (ix_img == px) then
      if (bc_right == BC_PERIODIC) then
         right_img = image_id(1, iy_img)
      else
         right_img = 0
      end if
   else
      right_img = image_id(ix_img + 1, iy_img)
   end if


   if (iy_img == 1) then
      if (bc_bottom == BC_PERIODIC) then
         down_img = image_id(ix_img, py)
      else
         down_img = 0
      end if
   else
      down_img = image_id(ix_img, iy_img - 1)
   end if


   if (iy_img == py) then
      if (bc_top == BC_PERIODIC) then
         up_img = image_id(ix_img, 1)
      else
         up_img = 0
      end if
   else
      up_img = image_id(ix_img, iy_img + 1)
   end if

    sync all

    if (my_img == 1) then
       print *, "LBM domain decomposition"
       print *, "Global grid: ", nx, " x ", ny
       print *, "Images:      ", px, " x ", py
       print *, "Local grid range: ", nx/px, "..", nx_alloc, " x ", &
         ny/py, "..", ny_alloc
    end if

  end subroutine setup_images


  subroutine choose_image_grid(n, px_out, py_out)
    integer, intent(in)  :: n
    integer, intent(out) :: px_out, py_out

    integer :: p, candidate_py
    real :: score, best_score

    ! Choose the factorization that produces the most square local blocks.
    px_out = 0
    py_out = 0
    best_score = huge(1.0)
    do p = 1, n
       if (mod(n, p) /= 0) cycle
       candidate_py = n/p
       if (p > nx .or. candidate_py > ny) cycle
       score = abs(real(nx)/real(p) - real(ny)/real(candidate_py))
       if (score < best_score) then
          best_score = score
          px_out = p
          py_out = candidate_py
       end if
    end do

    if (px_out == 0) then
       if (this_image() == 1) print *, &
         "Error: more images than can be assigned non-empty lattice blocks."
       error stop
    end if
  end subroutine choose_image_grid


  integer function block_size(global_size, parts, coordinate)
    integer, intent(in) :: global_size, parts, coordinate

    ! Return one block size from a quotient/remainder partition.
    block_size = global_size/parts
    if (coordinate <= mod(global_size,parts)) block_size = block_size + 1
  end function block_size


  integer function block_offset(global_size, parts, coordinate)
    integer, intent(in) :: global_size, parts, coordinate

    ! Return the zero-based global offset of one block.
    block_offset = (coordinate-1)*(global_size/parts) + &
      min(coordinate-1,mod(global_size,parts))
  end function block_offset


  integer function image_id(ix, iy)
    integer, intent(in) :: ix, iy

    ! Convert image-grid coordinates to a coarray image number.
    if (ix < 1 .or. ix > px .or. iy < 1 .or. iy > py) then
       image_id = 0
    else
       image_id = ix + (iy - 1) * px
    end if

  end function image_id



  subroutine communicate_ghost_cells(f)
  use lbm_params, only: dp, q
  implicit none

  real(dp), intent(inout) :: f(q, 0:nx_alloc+1, 0:ny_alloc+1)[*]
  integer :: remote_nx, remote_ny
  integer :: neighbour_images(4), neighbour_count

  ! Exchange post-collision edge cells for local pull streaming.
  neighbour_count = 0
  call add_neighbour(left_img,neighbour_images,neighbour_count)
  call add_neighbour(right_img,neighbour_images,neighbour_count)
  if (neighbour_count > 0) sync images(neighbour_images(:neighbour_count))

  if (left_img /= 0) then
     remote_nx = block_size(nx,px,merge(px,ix_img-1,ix_img == 1))
     f(:, 0, 1:ny_loc) = f(:, remote_nx, 1:ny_loc)[left_img]
  end if

  if (right_img /= 0) then
     f(:, nx_loc+1, 1:ny_loc) = f(:, 1, 1:ny_loc)[right_img]
  end if

  neighbour_count = 0
  call add_neighbour(left_img,neighbour_images,neighbour_count)
  call add_neighbour(right_img,neighbour_images,neighbour_count)
  call add_neighbour(down_img,neighbour_images,neighbour_count)
  call add_neighbour(up_img,neighbour_images,neighbour_count)
  if (neighbour_count > 0) sync images(neighbour_images(:neighbour_count))

  if (down_img /= 0) then
     remote_ny = block_size(ny,py,merge(py,iy_img-1,iy_img == 1))
     f(:, 0:nx_loc+1, 0) = f(:, 0:nx_loc+1, remote_ny)[down_img]
  end if

  if (up_img /= 0) then
     f(:, 0:nx_loc+1, ny_loc+1) = f(:, 0:nx_loc+1, 1)[up_img]
  end if

  neighbour_count = 0
  call add_neighbour(down_img,neighbour_images,neighbour_count)
  call add_neighbour(up_img,neighbour_images,neighbour_count)
  if (neighbour_count > 0) sync images(neighbour_images(:neighbour_count))

end subroutine communicate_ghost_cells


subroutine add_neighbour(image,images,count)
  implicit none

  integer, intent(in) :: image
  integer, intent(inout) :: images(:)
  integer, intent(inout) :: count

  ! Add a unique remote image to a synchronization list.
  if (image == 0 .or. image == my_img) return
  if (count > 0) then
     if (any(images(:count) == image)) return
  end if

  count = count + 1
  images(count) = image
end subroutine add_neighbour

end module domain_decomposition
