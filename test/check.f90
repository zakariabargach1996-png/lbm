program check
  use lbm_params, only: dp, q, cx, cy, opp, w, cs2, default_omega
  implicit none
  integer :: direction
  real(dp), parameter :: tolerance = 100.0_dp*epsilon(1.0_dp)

  ! Check the D2Q9 quadrature and opposite-direction identities.
  call assert_close(sum(w),1.0_dp,"D2Q9 weights sum to one")
  call assert_close(sum(w*real(cx,dp)),0.0_dp,"zero x momentum")
  call assert_close(sum(w*real(cy,dp)),0.0_dp,"zero y momentum")
  call assert_close(sum(w*real(cx*cx,dp)),cs2,"x second moment")
  call assert_close(sum(w*real(cy*cy,dp)),cs2,"y second moment")
  call assert_close(sum(w*real(cx*cy,dp)),0.0_dp,"mixed second moment")

  do direction = 1, q
    if (opp(opp(direction)) /= direction .or. &
        cx(opp(direction)) /= -cx(direction) .or. &
        cy(opp(direction)) /= -cy(direction)) then
      error stop "D2Q9 opposite-direction table is inconsistent"
    end if
  end do

  if (default_omega <= 0.0_dp .or. default_omega > 1.7_dp) then
    error stop "default omega lies outside the required stable range"
  end if
  print '(a)', "D2Q9 lattice and parameter tests passed."

contains

  subroutine assert_close(actual,expected,label)
    real(dp), intent(in) :: actual, expected
    character(len=*), intent(in) :: label
    ! Fail when two scalar values differ beyond the test tolerance.
    if (abs(actual-expected) > tolerance) then
      print '(a,2(1x,es24.16))', trim(label)//" failed:",actual,expected
      error stop 1
    end if
  end subroutine assert_close
end program check
