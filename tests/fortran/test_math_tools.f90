program test_math_tools
  use healpix_types
  use math_tools, only: mean, variance, tsum, corr_erf
  implicit none

  integer(i4b) :: failures
  real(dp), dimension(4) :: sample
  real(dp), dimension(3) :: x
  real(dp), dimension(3) :: y

  failures = 0
  sample = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
  x = [0.0_dp, 1.0_dp, 2.0_dp]
  y = [0.0_dp, 1.0_dp, 0.0_dp]

  call expect_close(mean(sample), 2.5_dp, 1.0e-12_dp, "mean", failures)
  call expect_close(variance(sample), 5.0_dp / 3.0_dp, 1.0e-12_dp, "variance", failures)
  call expect_close(tsum(x, y), 1.0_dp, 1.0e-12_dp, "tsum", failures)

  call expect_close(corr_erf(-1.0_dp), erf(-1.0_dp), 2.0e-7_dp, "corr_erf(-1)", failures)
  call expect_close(corr_erf(0.0_dp), erf(0.0_dp), 2.0e-7_dp, "corr_erf(0)", failures)
  call expect_close(corr_erf(1.0_dp), erf(1.0_dp), 2.0e-7_dp, "corr_erf(1)", failures)

  if (failures /= 0) then
    write(*, '(A,I0)') "test_math_tools failures: ", failures
    stop 1
  end if

  write(*, '(A)') "test_math_tools passed"

contains

  subroutine expect_close(actual, expected, tol, label, failures)
    real(dp), intent(in) :: actual, expected, tol
    character(len=*), intent(in) :: label
    integer(i4b), intent(inout) :: failures

    if (abs(actual - expected) > tol) then
      failures = failures + 1
      write(*, '(A,A,A,ES12.5,A,ES12.5)') "FAIL: ", trim(label), " expected=", expected, " actual=", actual
    end if
  end subroutine expect_close

end program test_math_tools
