module domain_decomposition
  use lbm_params
  implicit none

  ! This module decomposes the global nx-by-ny lattice over a two-dimensional
  ! grid of Fortran coarray images.  Image coordinates are one-based, as are
  ! physical local-cell indices.  The local storage convention is:
  !
  !     0        : lower/left ghost layer
  !     1:n*_loc : physical cells owned by this image
  !     n*_loc+1 : upper/right ghost layer
  !
  ! nx_alloc/ny_alloc are the largest local block sizes on any image.  Coarrays
  ! must have identical cobounds and array bounds on all images, so smaller
  ! blocks allocate the same extent but leave their extra interior slots unused.

  integer :: nimgs
  integer :: my_img

  integer :: px, py          ! Number of image blocks in global x and y.
  integer :: ix_img, iy_img  ! Coordinates of this image in that image grid.

  integer :: nx_loc, ny_loc     ! Owned physical cells; ghosts are excluded.
  integer :: nx_alloc, ny_alloc ! Common maximum block sizes for allocation.
  integer :: x_offset, y_offset ! Zero-based global origin of this block.

  ! A neighbour id of zero is a sentinel for a physical domain edge.  A
  ! periodic global edge instead points to the image on the opposite side.
  integer :: left_img, right_img
  integer :: down_img, up_img

contains

  subroutine setup_images()
    implicit none

    ! Coarray images are numbered 1..num_images().  choose_image_grid factors
    ! that one-dimensional set into a px-by-py Cartesian topology.
    nimgs  = num_images()
    my_img = this_image()

    call choose_image_grid(nimgs, px, py)

    ! Row-major mapping: x changes fastest in the coarray image number.
    ! image_id below implements the inverse transformation.
    ix_img = mod(my_img - 1, px) + 1
    iy_img = (my_img - 1) / px + 1

    ! If nx or ny is not divisible by the number of blocks, the first remainder
    ! blocks receive one extra cell.  Offsets account for those larger blocks.
    nx_loc = block_size(nx, px, ix_img)
    ny_loc = block_size(ny, py, iy_img)
    nx_alloc = (nx + px - 1) / px
    ny_alloc = (ny + py - 1) / py
    x_offset = block_offset(nx, px, ix_img)
    y_offset = block_offset(ny, py, iy_img)



   ! Build the horizontal neighbour map.  At a periodic global edge, wrap to
   ! the last/first image in the same image row.  At a wall, use zero so halo
   ! exchange will leave that ghost layer untouched for boundary reconstruction.
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


   ! Build the corresponding vertical neighbour map.  Internal block edges
   ! always have a real neighbour regardless of the physical boundary types.
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

    ! All images must finish publishing the same decomposition state before
    ! the simulation allocates coarrays and starts remote accesses.
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

    ! Consider every factorisation n = px*py that gives every image at least
    ! one lattice cell.  A nearly square local block has a smaller perimeter,
    ! which generally reduces the amount of halo data relative to useful work.
    px_out = 0
    py_out = 0
    best_score = huge(1.0)
    do p = 1, n
       if (mod(n, p) /= 0) cycle
       candidate_py = n/p
       if (p > nx .or. candidate_py > ny) cycle
       ! Prefer nearly square local domains, reducing halo surface area.
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

    ! Quotient/remainder partitioning keeps block sizes within one cell of one
    ! another.  The lower-coordinate blocks own the remainder cells.
    block_size = global_size/parts
    if (coordinate <= mod(global_size,parts)) block_size = block_size + 1
  end function block_size


  integer function block_offset(global_size, parts, coordinate)
    integer, intent(in) :: global_size, parts, coordinate

    ! Number of cells owned by all lower-coordinate blocks: their base-size
    ! cells plus one extra cell for each preceding remainder block.
    block_offset = (coordinate-1)*(global_size/parts) + &
      min(coordinate-1,mod(global_size,parts))
  end function block_offset


  integer function image_id(ix, iy)
    integer, intent(in) :: ix, iy

    ! Convert Cartesian image coordinates back to a coarray image number.
    ! Coordinates outside the image grid return the physical-boundary sentinel.
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

  ! f is a coarray: each image can read another image's local f by appending a
  ! coindex, for example f(...)[left_img].  Only f needs to be coarray storage;
  ! macroscopic fields and the next-time-step buffer are strictly local.
  !
  ! The solver uses pull streaming:
  !   f_next(q,i,j) = f(q,i-cx(q),j-cy(q)).
  ! Therefore every physical cell can be streamed locally once the one-cell
  ! halo holds the neighbouring post-collision populations.  One layer is
  ! sufficient because every D2Q9 velocity component is -1, 0, or +1.
  !
  ! The initial barrier ensures every image has completed collision before any
  ! image reads a neighbour's post-collision boundary.  This is deliberately a
  ! simple, correctness-first synchronization scheme.
  sync all

  ! Fill the left and right ghost columns from the adjacent images' physical
  ! edge columns.  A zero neighbour means a solid physical edge; its ghost
  ! values are not used as boundary data and are corrected after streaming.

  if (left_img /= 0) then
     ! The left neighbour may own a different number of x cells.  Its physical
     ! right edge is remote_nx, not necessarily this image's nx_loc.
     remote_nx = block_size(nx,px,merge(px,ix_img-1,ix_img == 1))
     f(:, 0, 1:ny_loc) = f(:, remote_nx, 1:ny_loc)[left_img]
  end if

  if (right_img /= 0) then
     f(:, nx_loc+1, 1:ny_loc) = f(:, 1, 1:ny_loc)[right_img]
  end if

  ! Complete all horizontal reads before exchanging rows.  The row exchange
  ! includes x=0 and x=nx_loc+1, so it also propagates the freshly received
  ! horizontal ghost values into diagonal corner ghosts.  Those corners supply
  ! diagonal D2Q9 populations when both x and y cross an image boundary.
  sync all

  ! Fill bottom and top ghost rows.  As above, the lower neighbour's physical
  ! top row depends on that neighbour's actual (possibly uneven) block size.

  if (down_img /= 0) then
     remote_ny = block_size(ny,py,merge(py,iy_img-1,iy_img == 1))
     f(:, 0:nx_loc+1, 0) = f(:, 0:nx_loc+1, remote_ny)[down_img]
  end if

  if (up_img /= 0) then
     f(:, 0:nx_loc+1, ny_loc+1) = f(:, 0:nx_loc+1, 1)[up_img]
  end if

  ! No image may start overwriting f in the next collision step while another
  ! image could still be reading it remotely for this exchange.
  sync all

end subroutine communicate_ghost_cells

end module domain_decomposition
