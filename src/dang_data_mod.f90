module dang_data_mod
  use healpix_types
  use utility_mod
  implicit none

  type, public                                :: data

     ! Storage
    character(len=255)                        :: label
    
    integer(i4b)                              :: npix
    integer(i4b)                              :: nfgs
    integer(i4b)                              :: nmaps

    real(dp), allocatable, dimension(:,:,:)   :: sig_map
    real(dp), allocatable, dimension(:,:,:)   :: rms_map
    real(dp), allocatable, dimension(:,:,:)   :: res_map
    real(dp), allocatable, dimension(:,:)     :: chi_map
    real(dp), allocatable, dimension(:,:,:,:) :: fg_map

    real(dp), allocatable, dimension(:,:,:)   :: temps

    real(dp), allocatable, dimension(:,:,:)   :: temp_amps
  end type data

contains

  subroutine init_fg_map(self,npix,nmaps,nbands,nfgs)
    implicit none
    type(data)               :: self
    integer(i4b), intent(in) :: npix, nmaps, nbands, nfgs
    
    allocate(self%fg_map(0:npix-1,nmaps,nbands,nfgs))
    
  end subroutine init_fg_map

  subroutine init_data_maps(self,npix,nmaps,nbands)
    implicit none
    type(data)               :: self
    integer(i4b), intent(in) :: npix, nmaps, nbands

    allocate(self%sig_map(0:npix-1,nmaps,nbands))
    allocate(self%rms_map(0:npix-1,nmaps,nbands))
    allocate(self%res_map(0:npix-1,nmaps,nbands))
    allocate(self%chi_map(0:npix-1,nmaps))

  end subroutine init_data_maps

  subroutine init_template(self,npix,nmaps,ntemp)
    implicit none
    type(data)               :: self
    integer(i4b), intent(in) :: npix, nmaps, ntemp

    allocate(self%temps(0:npix-1,nmaps,ntemp))

  end subroutine init_template

  subroutine init_temp_amps(self,nbands,nmaps,ntemp)
    implicit none
    type(data)               :: self
    integer(i4b), intent(in) :: nbands, nmaps, ntemp

    allocate(self%temp_amps(nbands,nmaps,ntemp))

  end subroutine init_temp_amps


end module dang_data_mod
