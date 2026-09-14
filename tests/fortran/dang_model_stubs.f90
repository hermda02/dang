module dang_param_mod
  use healpix_types
  implicit none

  type dang_params
     integer(i4b) :: dummy = 0
  end type dang_params
end module dang_param_mod

module dang_data_mod
  use healpix_types
  implicit none

  type dang_data
     integer(i4b) :: npix = 0
     real(dp), allocatable, dimension(:,:,:) :: rms_map
  end type dang_data
end module dang_data_mod

module dang_bp_mod
  use healpix_types
  implicit none

  type, abstract :: bandinfo
     integer(i4b) :: n = 0
     real(dp), allocatable, dimension(:) :: nu0, tau0
   contains
     procedure(integrate_iface), deferred :: integrate
  end type bandinfo

  abstract interface
     function integrate_iface(self, signal) result(integrated_signal)
       import :: bandinfo, dp
       class(bandinfo), intent(in) :: self
       real(dp), dimension(:), intent(in) :: signal
       real(dp) :: integrated_signal
     end function integrate_iface
  end interface
end module dang_bp_mod

module dang_component_mod
  use healpix_types
  implicit none

  type, abstract :: dang_comp
     real(dp) :: nu_ref = 0.0_dp
     real(dp), allocatable, dimension(:,:) :: uni_prior
   contains
     procedure(evalSED), deferred :: S
  end type dang_comp

  abstract interface
     function evalSED(self, nu, band, pol, theta, pixel) result(eval_sed)
       import :: dang_comp, dp, i4b
       class(dang_comp), intent(in) :: self
       real(dp),                intent(in), optional :: nu
       integer(i4b),            intent(in), optional :: band
       integer(i4b),            intent(in), optional :: pixel
       integer(i4b),            intent(in)           :: pol
       real(dp), dimension(1:), intent(in), optional :: theta
       real(dp)                                      :: eval_sed
     end function evalSED
  end interface
end module dang_component_mod
