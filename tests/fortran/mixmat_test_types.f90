module mixmat_test_types
  use healpix_types
  use dang_component_mod, only: dang_comp
  use dang_bp_mod, only: bandinfo
  implicit none

  type, extends(bandinfo) :: test_band
   contains
    procedure :: integrate => test_integrate
  end type test_band

  type, extends(dang_comp) :: test_comp
   contains
    procedure :: S => test_s
  end type test_comp

contains

  function test_integrate(self, signal) result(integrated_signal)
    class(test_band), intent(in) :: self
    real(dp), dimension(:), intent(in) :: signal
    real(dp) :: integrated_signal

    integrated_signal = sum(self%tau0 * signal)
  end function test_integrate

  function test_s(self, nu, band, pol, theta, pixel) result(eval_sed)
    class(test_comp), intent(in) :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pixel
    integer(i4b),            intent(in)           :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: eval_sed

    eval_sed = 0.0_dp
    if (present(band)) eval_sed = eval_sed + 0.0_dp * real(band, dp)
    if (present(pixel)) eval_sed = eval_sed + 0.0_dp * real(pixel, dp)
    eval_sed = eval_sed + 0.0_dp * real(pol, dp) + 0.0_dp * self%nu_ref

    if (.not. present(theta) .or. .not. present(nu)) then
      return
    else if (size(theta) == 1) then
      eval_sed = theta(1) * nu + 2.0_dp
    else
      eval_sed = theta(1) + 2.0_dp * theta(2) + 0.5_dp * nu
    end if
  end function test_s

end module mixmat_test_types
