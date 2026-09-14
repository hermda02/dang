program test_dang_mixmat_mod
  use healpix_types
  use dang_mixmat_1d_mod, only: dang_mixmat_1d
  use dang_mixmat_2d_mod, only: dang_mixmat_2d
  use mixmat_test_types, only: test_band, test_comp
  implicit none

  integer(i4b) :: failures
  type(test_band), target :: bp
  type(test_comp), target :: comp
  class(dang_mixmat_1d), pointer :: m1
  class(dang_mixmat_2d), pointer :: m2
  real(dp), dimension(1) :: theta1
  real(dp), dimension(2) :: theta2
  real(dp) :: nu_bar

  failures = 0

  call init_inputs(comp, bp, nu_bar)

  m1 => dang_mixmat_1d(comp, bp, 1_i4b)
  theta1 = [0.9_dp]
  call expect_close(m1%eval(theta1), theta1(1) * nu_bar + 2.0_dp, 2.0e-3_dp, "mixmat_1d integral", failures)
  call expect_close(m1%eval_dI(theta1, 1_i4b), nu_bar, 2.0e-3_dp, "mixmat_1d derivative", failures)

  m2 => dang_mixmat_2d(comp, bp, 1_i4b)
  theta2 = [0.4_dp, -0.6_dp]
  call expect_close(m2%eval(theta2), theta2(1) + 2.0_dp * theta2(2) + 0.5_dp * nu_bar, 1.0e-2_dp, "mixmat_2d integral", failures)
  call expect_close(m2%eval_dI(theta2, 1_i4b), 1.0_dp, 5.0e-3_dp, "mixmat_2d d/dtheta1", failures)
  call expect_close(m2%eval_dI(theta2, 2_i4b), 2.0_dp, 5.0e-3_dp, "mixmat_2d d/dtheta2", failures)

  if (failures /= 0) then
    write(*, '(A,I0)') "test_dang_mixmat_mod failures: ", failures
    stop 1
  end if

  write(*, '(A)') "test_dang_mixmat_mod passed"

contains

  subroutine init_inputs(comp, bp, nu_bar)
    type(test_comp), intent(inout) :: comp
    type(test_band), intent(inout) :: bp
    real(dp), intent(out) :: nu_bar

    allocate(comp%uni_prior(2, 2))
    comp%uni_prior(1, :) = [0.2_dp, 1.4_dp]
    comp%uni_prior(2, :) = [-1.0_dp, 1.0_dp]
    comp%nu_ref = 10.0_dp

    bp%n = 3
    allocate(bp%nu0(3), bp%tau0(3))
    bp%nu0 = [1.0_dp, 2.0_dp, 4.0_dp]
    bp%tau0 = [0.2_dp, 0.3_dp, 0.5_dp]
    nu_bar = sum(bp%tau0 * bp%nu0)
  end subroutine init_inputs

  subroutine expect_close(actual, expected, tol, label, failures)
    real(dp), intent(in) :: actual, expected, tol
    character(len=*), intent(in) :: label
    integer(i4b), intent(inout) :: failures

    if (abs(actual - expected) > tol) then
      failures = failures + 1
      write(*, '(A,A,A,ES12.5,A,ES12.5)') "FAIL: ", trim(label), " expected=", expected, " actual=", actual
    end if
  end subroutine expect_close

end program test_dang_mixmat_mod
