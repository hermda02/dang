program test_locate_mod
  use healpix_types
  use locate_mod, only: locate
  implicit none

  integer(i4b) :: failures
  integer(i4b), dimension(5) :: xi
  real(dp), dimension(5) :: xr

  failures = 0
  xi = [1_i4b, 3_i4b, 5_i4b, 7_i4b, 9_i4b]
  xr = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 16.0_dp]

  call expect_equal_i4(locate(xi, 1_i4b), 1_i4b, "locate_int lower edge", failures)
  call expect_equal_i4(locate(xi, 9_i4b), 4_i4b, "locate_int upper edge", failures)
  call expect_equal_i4(locate(xi, 6_i4b), 3_i4b, "locate_int midpoint", failures)

  call expect_equal_i4(locate(xr, 1.0_dp), 1_i4b, "locate_dp lower edge", failures)
  call expect_equal_i4(locate(xr, 16.0_dp), 4_i4b, "locate_dp upper edge", failures)
  call expect_equal_i4(locate(xr, 3.0_dp), 2_i4b, "locate_dp midpoint", failures)

  if (failures /= 0) then
    write(*, '(A,I0)') "test_locate_mod failures: ", failures
    stop 1
  end if

  write(*, '(A)') "test_locate_mod passed"

contains

  subroutine expect_equal_i4(actual, expected, label, failures)
    integer(i4b), intent(in) :: actual, expected
    character(len=*), intent(in) :: label
    integer(i4b), intent(inout) :: failures

    if (actual /= expected) then
      failures = failures + 1
      write(*, '(A,A,A,I0,A,I0)') "FAIL: ", trim(label), " expected=", expected, " actual=", actual
    end if
  end subroutine expect_equal_i4

end program test_locate_mod
