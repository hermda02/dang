module foreground_mod
  use healpix_types
  use init_mod
  implicit none
  !// Make everything not specified as public invisible from outside the
  !// module
  private
  type, public          :: fg_comp
    !// Internal variables for this type
    real(dp)            :: nu_ref      ! Foreground reference frequency
    integer             :: loc         ! Reference frequency band (nearest)
    real(dp)            :: p(2)        ! Foreground Parameters
    real(dp)            :: gauss_mean  ! Parameter sampling mean
    real(dp)            :: gauss_std   ! Parameter sampling standard dev
    real(dp)            :: uni_min     ! Uniform minimum (won't accept samples below this value)
    real(dp)            :: uni_max     ! Uniform maximum (won't accept samples above this value)
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

  end module foreground_mod