program test_dang_util_mod
  use healpix_types
  use dang_util_mod, only: eval_normal_prior, tsum, count_delimits, delimit_string
  use dang_util_mod, only: return_poltype_flag, mask_avg, mask_sum, npix
  implicit none

  integer(i4b) :: failures

  failures = 0

  call test_eval_normal_prior_formula(failures)
  call test_tsum_linear_exact(failures)
  call test_delimit_helpers(failures)
  call test_return_poltype_flag(failures)
  call test_mask_aggregates(failures)

  if (failures /= 0) then
    write(*, '(A,I0)') "test_dang_util_mod failures: ", failures
    stop 1
  end if

  write(*, '(A)') "test_dang_util_mod passed"

contains

  subroutine test_eval_normal_prior_formula(failures)
    integer(i4b), intent(inout) :: failures
    real(dp) :: p, expected

    p = eval_normal_prior(1.5_dp, 1.0_dp, 2.0_dp)
    expected = exp(-((1.5_dp - 1.0_dp)**2) / (2.0_dp * 4.0_dp)) / (2.0_dp * sqrt(2.0_dp * pi))
    call expect_close(p, expected, 1.0e-14_dp, "eval_normal_prior formula", failures)
  end subroutine test_eval_normal_prior_formula

  subroutine test_tsum_linear_exact(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(5) :: x, y
    real(dp) :: exact

    x = [0.0_dp, 0.5_dp, 1.5_dp, 2.0_dp, 3.0_dp]
    y = 3.0_dp * x + 2.0_dp
    exact = 0.5_dp * 3.0_dp * 3.0_dp**2 + 2.0_dp * 3.0_dp
    call expect_close(tsum(x, y), exact, 1.0e-12_dp, "tsum linear exactness", failures)
  end subroutine test_tsum_linear_exact

  subroutine test_delimit_helpers(failures)
    integer(i4b), intent(inout) :: failures
    character(len=5), dimension(3) :: tokens

    if (count_delimits("T,Q+U,U", ",") /= 2_i4b) then
      failures = failures + 1
      write(*, '(A)') "FAIL: count_delimits"
    end if

    call delimit_string("T,Q+U,U", ",", tokens)
    if (tokens(1) /= "T" .or. tokens(2) /= "Q+U" .or. tokens(3) /= "U") then
      failures = failures + 1
      write(*, '(A)') "FAIL: delimit_string"
    end if
  end subroutine test_delimit_helpers

  subroutine test_return_poltype_flag(failures)
    integer(i4b), intent(inout) :: failures
    integer(i4b), allocatable, dimension(:) :: flag

    flag = return_poltype_flag("T,Q+U     ")
    if (size(flag) /= 2) then
      failures = failures + 1
      write(*, '(A)') "FAIL: return_poltype_flag size"
      return
    end if
    if (flag(1) /= 1_i4b .or. flag(2) /= 8_i4b) then
      failures = failures + 1
      write(*, '(A)') "FAIL: return_poltype_flag values"
    end if
  end subroutine test_return_poltype_flag

  subroutine test_mask_aggregates(failures)
    integer(i4b), intent(inout) :: failures
    real(dp), dimension(4) :: arr, mask
    real(dp), parameter :: missing = -1.6375d30

    npix = 4
    arr = [10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp]
    mask = [1.0_dp, 0.0_dp, missing, 1.0_dp]

    call expect_close(mask_sum(arr, mask), 50.0_dp, 1.0e-12_dp, "mask_sum ignores zero/miss", failures)
    call expect_close(mask_avg(arr, mask), 25.0_dp, 1.0e-12_dp, "mask_avg ignores zero/miss", failures)
  end subroutine test_mask_aggregates

  subroutine expect_close(actual, expected, tol, label, failures)
    real(dp), intent(in) :: actual, expected, tol
    character(len=*), intent(in) :: label
    integer(i4b), intent(inout) :: failures
    if (abs(actual - expected) > tol) then
      failures = failures + 1
      write(*, '(A,A,A,ES12.5,A,ES12.5)') "FAIL: ", trim(label), " expected=", expected, " actual=", actual
    end if
  end subroutine expect_close

end program test_dang_util_mod
