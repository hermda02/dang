program test_spline_2d_mod
  use healpix_types
  use spline_2D_mod, only: splie2_full_precomp, splin2_full_precomp, splie2, splin2
  implicit none

  integer(i4b) :: failures

  failures = 0

  call test_affine_exactness(failures)
  call test_grid_node_reproduction(failures)
  call test_precomp_matches_direct(failures)

  if (failures /= 0) then
    write(*, '(A,I0)') "test_spline_2d_mod failures: ", failures
    stop 1
  end if

  write(*, '(A)') "test_spline_2d_mod passed"

contains

  subroutine test_affine_exactness(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(4) :: x, y
    real(dp), dimension(4, 4) :: f, y2
    real(dp), dimension(4, 4, 4, 4) :: coeff
    real(dp), dimension(4) :: xs, ys
    integer(i4b) :: i, j

    x = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
    y = [0.0_dp, 1.5_dp, 3.0_dp, 4.5_dp]

    do i = 1, size(x)
      do j = 1, size(y)
        f(i, j) = affine_f(x(i), y(j))
      end do
    end do

    call splie2_full_precomp(x, y, f, coeff)

    xs = [0.3_dp, 1.1_dp, 2.2_dp, 2.7_dp]
    ys = [0.9_dp, 2.1_dp, 3.3_dp, 4.2_dp]
    do i = 1, size(xs)
      call expect_close(splin2_full_precomp(x, y, coeff, xs(i), ys(i)), affine_f(xs(i), ys(i)), 1.0e-10_dp, &
        "full-precomp affine exact", failures)
    end do

    call splie2(x, y, f, y2)
    do i = 1, size(xs)
      call expect_close(splin2(x, y, f, y2, xs(i), ys(i)), affine_f(xs(i), ys(i)), 1.0e-10_dp, &
        "direct affine exact", failures)
    end do
  end subroutine test_affine_exactness

  subroutine test_grid_node_reproduction(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(4) :: x, y
    real(dp), dimension(4, 4) :: f, y2
    real(dp), dimension(4, 4, 4, 4) :: coeff
    integer(i4b) :: i, j

    x = [0.0_dp, 0.7_dp, 1.7_dp, 2.5_dp]
    y = [0.0_dp, 0.8_dp, 1.9_dp, 3.1_dp]

    do i = 1, size(x)
      do j = 1, size(y)
        f(i, j) = smooth_f(x(i), y(j))
      end do
    end do

    call splie2_full_precomp(x, y, f, coeff)
    call splie2(x, y, f, y2)

    do i = 1, size(x)
      do j = 1, size(y)
        call expect_close(splin2_full_precomp(x, y, coeff, x(i), y(j)), f(i, j), 1.0e-9_dp, "full-precomp nodes", failures)
        call expect_close(splin2(x, y, f, y2, x(i), y(j)), f(i, j), 1.0e-9_dp, "direct nodes", failures)
      end do
    end do
  end subroutine test_grid_node_reproduction

  subroutine test_precomp_matches_direct(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(5) :: x, y
    real(dp), dimension(5, 5) :: f, y2
    real(dp), dimension(4, 4, 5, 5) :: coeff
    real(dp), dimension(3) :: xs, ys
    integer(i4b) :: i, j
    real(dp) :: fp, fd

    x = [0.0_dp, 0.6_dp, 1.2_dp, 2.0_dp, 2.8_dp]
    y = [0.0_dp, 0.5_dp, 1.4_dp, 2.4_dp, 3.0_dp]

    do i = 1, size(x)
      do j = 1, size(y)
        f(i, j) = smooth_f(x(i), y(j))
      end do
    end do

    call splie2_full_precomp(x, y, f, coeff)
    call splie2(x, y, f, y2)

    xs = [0.45_dp, 1.35_dp, 2.45_dp]
    ys = [0.35_dp, 1.75_dp, 2.65_dp]
    do i = 1, size(xs)
      fp = splin2_full_precomp(x, y, coeff, xs(i), ys(i))
      fd = splin2(x, y, f, y2, xs(i), ys(i))
      call expect_close(fp, fd, 4.0e-2_dp, "full-precomp vs direct", failures)
    end do
  end subroutine test_precomp_matches_direct

  pure real(dp) function affine_f(x, y)
    real(dp), intent(in) :: x, y
    affine_f = 1.25_dp + 2.0_dp * x - 0.5_dp * y
  end function affine_f

  pure real(dp) function smooth_f(x, y)
    real(dp), intent(in) :: x, y
    smooth_f = x**2 + 0.3_dp * x * y + sin(y)
  end function smooth_f

  subroutine expect_close(actual, expected, tol, label, failures)
    real(dp), intent(in) :: actual, expected, tol
    character(len=*), intent(in) :: label
    integer(i4b), intent(inout) :: failures

    if (abs(actual - expected) > tol) then
      failures = failures + 1
      write(*, '(A,A,A,ES12.5,A,ES12.5)') "FAIL: ", trim(label), " expected=", expected, " actual=", actual
    end if
  end subroutine expect_close

end program test_spline_2d_mod
