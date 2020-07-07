module foreground_mod
  use healpix_types
  use init_mod
  implicit none
  !// Make everything not specified as public invisible from outside the
  !// module
  private
  type, public          :: fg_comp
    !// Internal variables for this type
    real(dp)            :: nu_ref        ! Foreground reference frequency
    integer             :: loc           ! Reference frequency band (nearest)
    real(dp)            :: p(2)          ! Foreground Parameters
    real(dp)            :: gauss_mean(2) ! Parameter sampling mean
    real(dp)            :: gauss_std(2)  ! Parameter sampling standard dev
    real(dp)            :: uni_min(2)    ! Uniform minimum (won't accept samples below this value)
    real(dp)            :: uni_max(2)    ! Uniform maximum (won't accept samples above this value)
    character(len=255)  :: type
  contains
    !// the procedures to load data, to write data to the screen,
    !// to clear the contents of an object and to copy data from one object
    !// to another
    procedure           :: get_params => get_comp_params
    procedure           :: get_nu_ref => get_comp_nu_ref
    procedure           :: get_type => get_comp_type
    procedure           :: get_gauss_mean => get_comp_gauss_mean
    procedure           :: get_gauss_std  => get_comp_gauss_std
  end type fg_comp

contains
  subroutine get_comp_nu_ref(self, freq)
    class(fg_comp)            :: self
    real(dp), intent(in)      :: freq

    self%nu_ref = freq
  end subroutine get_comp_nu_ref

  subroutine get_comp_type(self, typ)
    class(fg_comp)               :: self
    character(len=*), intent(in) :: typ

    self%type = typ
  end subroutine get_comp_type

  subroutine get_comp_params(self,params)
    class(fg_comp)               :: self
    real(dp), intent(in)         :: params(:)

    self%p = params
  end subroutine get_comp_params

  subroutine get_comp_gauss_mean(self,gmean)
    class(fg_comp)             :: self
    real(dp), intent(in)       :: gmean

    self%gauss_mean = gmean
  end subroutine get_comp_gauss_mean

  subroutine get_comp_gauss_std(self,gstd)
    class(fg_comp)             :: self
    real(dp), intent(in)       :: gstd

    self%gauss_std = gstd
  end subroutine get_comp_gauss_std

  ! function compute_spectrum(self, freq) result(spectrum)
  !   implicit none
  !   class(fg_comp)                 :: self
  !   real(dp),           intent(in) :: freq
  !   real(dp)                       :: y, rj_cmb ! Conversion factor from K_{RJ} -> K_{CMB}
  !   real(dp)                       :: z!, compute_spectrum
  !   real(dp)                       :: spectrum

  !   y = h*(freq*1.0d9) / (k_B*T_CMB)
  !   rj_cmb = ((exp(y)-1)**2.d0)/(y**2.d0*exp(y))

  !   if (trim(self%type) == 'synch') then
  !     spectrum = (freq/self%nu_ref)**self%p(1) !*rj_cmb
  !   else if (trim(self%type) == 'dust') then
  !     z = h / (k_B*self%p(2))
  !     spectrum = (exp(z*self%nu_ref*1d9)-1.d0) / (exp(z*freq*1d9)-1.d0) * (freq/self%nu_ref)**(self%p(1)+1.d0)!*rj_cmb
  !   end if

  ! end function compute_spectrum

end module foreground_mod