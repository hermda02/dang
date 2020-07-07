module data_mod
  use healpix_types
  use init_mod
  implicit none
  !// Make everything not specified as public invisible from outside the
  !// module
  private
  !// Declare a  type called matrix
  type, public                              :: band
    !// Internal variables for this type
    real(dp)                                :: nu
    character(len=255)                      :: sig_map
    character(len=255)                      :: rms_map
    character(len=255)                      :: label
  contains
    !// the procedures to load data, to write data to the screen,
    !// to clear the contents of an object and to copy data from one object
    !// to another
    procedure                               :: get_nu => get_band_nu
    procedure                               :: get_map => get_band_file
    procedure                               :: get_rms => get_band_rms
    procedure                               :: get_label => get_band_label

  end type band

contains

    subroutine get_band_nu(self, freq)
        class(band)            :: self
        real(dp), intent(in)   :: freq

        self%nu = freq
    end subroutine get_band_nu

    subroutine get_band_file(self, filename)
      class(band)                  :: self
      character(len=*), intent(in) :: filename

      self%sig_map = filename

    end subroutine get_band_file

    subroutine get_band_rms(self, filename)
      class(band)                  :: self
      character(len=*), intent(in) :: filename

      self%rms_map = filename

    end subroutine get_band_rms

    subroutine get_band_label(self, name)
      class(band)                  :: self
      character(len=*), intent(in) :: name

      self%label = name

    end subroutine get_band_label

end module data_mod
