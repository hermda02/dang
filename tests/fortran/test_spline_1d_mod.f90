program test_spline_1d_mod
  use healpix_types
  use spline_1D_mod, only: spline_type, spline, splint, splint_multi, spline_plain, splint_plain, free_spline
  implicit none

  integer(i4b) :: failures
  type(spline_type) :: s
  real(dp), dimension(3) :: x_nodes, y_nodes, y2_nodes
  real(dp), dimension(3) :: xs, ys

  failures = 0

  x_nodes = [0.0_dp, 1.0_dp, 2.0_dp]
  y_nodes = 2.0_dp * x_nodes + 1.0_dp

  call spline_plain(x_nodes, y_nodes, 1.0d30, 1.0d30, y2_nodes)
  call expect_close_scalar(y2_nodes(1), 0.0_dp, 1.0e-12_dp, "spline_plain y2(1)", failures)
  call expect_close_scalar(y2_nodes(2), 0.0_dp, 1.0e-12_dp, "spline_plain y2(2)", failures)
  call expect_close_scalar(y2_nodes(3), 0.0_dp, 1.0e-12_dp, "spline_plain y2(3)", failures)
  call expect_close_scalar(splint_plain(x_nodes, y_nodes, y2_nodes, 0.5_dp), 2.0_dp, 1.0e-12_dp, "splint_plain x=0.5", failures)
  call expect_close_scalar(splint_plain(x_nodes, y_nodes, y2_nodes, 1.5_dp), 4.0_dp, 1.0e-12_dp, "splint_plain x=1.5", failures)

  call spline(s, x_nodes, y_nodes, regular=.true.)
  call expect_close_scalar(splint(s, 0.25_dp), 1.5_dp, 1.0e-12_dp, "splint simple regular", failures)
  call splint_multi(s, [0.25_dp, 0.75_dp, 1.75_dp], ys)
  call expect_close_scalar(ys(1), 1.5_dp, 1.0e-12_dp, "splint_multi 1", failures)
  call expect_close_scalar(ys(2), 2.5_dp, 1.0e-12_dp, "splint_multi 2", failures)
  call expect_close_scalar(ys(3), 4.5_dp, 1.0e-12_dp, "splint_multi 3", failures)

  call free_spline(s)
  xs = [0.0_dp, 2.0_dp, 5.0_dp]
  ys = [0.0_dp, 4.0_dp, 10.0_dp]
  call spline(s, xs, ys, linear=.true.)
  call expect_close_scalar(splint(s, 3.0_dp), 6.0_dp, 1.0e-12_dp, "splint simple linear", failures)
  call free_spline(s)

  if (failures /= 0) then
    write(*, '(A,I0)') "test_spline_1d_mod failures: ", failures
    stop 1
  end if

  write(*, '(A)') "test_spline_1d_mod passed"

contains

  subroutine expect_close_scalar(actual, expected, tol, label, failures)
    real(dp), intent(in) :: actual, expected, tol
    character(len=*), intent(in) :: label
    integer(i4b), intent(inout) :: failures

    if (abs(actual - expected) > tol) then
      failures = failures + 1
      write(*, '(A,A,A,ES12.5,A,ES12.5)') "FAIL: ", trim(label), " expected=", expected, " actual=", actual
    end if
  end subroutine expect_close_scalar

end program test_spline_1d_mod
