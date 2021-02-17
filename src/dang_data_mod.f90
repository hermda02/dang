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

    real(dp), allocatable, dimension(:)       :: gain
    real(dp), allocatable, dimension(:)       :: offset
    
    real(dp), allocatable, dimension(:,:)     :: masks

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

  subroutine init_mask_maps(self,npix,nmaps)
    implicit none
    type(data)               :: self
    integer(i4b), intent(in) :: npix, nmaps
    
    allocate(self%masks(0:npix-1,nmaps))
    
  end subroutine init_mask_maps

  subroutine init_data_maps(self,npix,nmaps,nbands)
    implicit none
    type(data)               :: self
    integer(i4b), intent(in) :: npix, nmaps, nbands

    allocate(self%sig_map(0:npix-1,nmaps,nbands))
    allocate(self%rms_map(0:npix-1,nmaps,nbands))
    allocate(self%res_map(0:npix-1,nmaps,nbands))
    allocate(self%chi_map(0:npix-1,nmaps))

    allocate(self%gain(nbands))
    allocate(self%offset(nbands))

    self%gain   = 1.d0
    self%offset = 0.d0

  end subroutine init_data_maps

end module dang_data_mod
