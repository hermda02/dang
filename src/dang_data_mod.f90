module dang_data_mod
  use healpix_types
  use dang_util_mod
  use dang_param_mod
  implicit none

  type, public                                :: data

     ! Storage
    character(len=255)                        :: label
    
    integer(i4b)                              :: npix
    integer(i4b)                              :: nfgs
    integer(i4b)                              :: nmaps

    real(dp), allocatable, dimension(:,:,:)   :: sig_map    ! data maps 
    real(dp), allocatable, dimension(:,:,:)   :: rms_map    ! noise maps
    real(dp), allocatable, dimension(:,:,:)   :: res_map    ! Residual maps
    real(dp), allocatable, dimension(:,:)     :: chi_map    ! Chisq map (for outputs)
    real(dp), allocatable, dimension(:,:,:,:) :: fg_map     ! Component maps

    real(dp), allocatable, dimension(:)       :: gain       ! Where band gains are stored
    real(dp), allocatable, dimension(:)       :: offset     ! Where band offsets are stored
    
    real(dp), allocatable, dimension(:,:)     :: masks      ! Where masks are stored

    real(dp), allocatable, dimension(:,:,:)   :: temps      ! Where template maps are stored   
    real(dp), allocatable, dimension(:,:,:)   :: temp_amps  ! Where template amplitudes are stored
    real(dp), allocatable, dimension(:,:)     :: temp_norm  ! Where we store template normalizations

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

  function compute_bnu_prime_RJ(nu)
    ! Assume that nu is in GHz
    implicit none
    real(dp), intent(in) :: nu
    real(dp)             :: compute_bnu_prime_RJ

    compute_bnu_prime_RJ = 2.d0*k_B*nu**2.d0/c**2.d0

  end function compute_bnu_prime_RJ

  function compute_bnu_prime(nu)
    ! Assume that nu is in GHz
    implicit none
    real(dp), intent(in) :: nu
    real(dp)             :: compute_bnu_prime, y

    y = h*nu/(k_B*T_CMB)
    compute_bnu_prime = (2.d0*h*nu**3)/(c**2.d0*(exp(y)-1))*(exp(y)/(exp(y)-1))*h*nu/(k_B*T_CMB**2)

  end function compute_bnu_prime

  function a2f(nu)
    ! [MJy/sr / uK_RJ]
    ! Assume that nu is in GHz
    implicit none
    real(dp), intent(in) :: nu
    real(dp)             :: a2f, y

    a2f = compute_bnu_prime_RJ(nu)

  end function a2f

  function a2t(nu)
    ! [uK_cmb/uK_RJ]
    ! Assume that nu is in GHz
    implicit none
    real(dp), intent(in) :: nu
    real(dp)             :: a2t, y

    y = (h*nu)/(k_B*T_CMB)

    a2t = (exp(y)-1.d0)**2/(y**2*exp(y))

  end function a2t

  function f2t(nu)
    ! [uK_cmb/MJysr-2]
    ! Assume that nu is in GHz
    implicit none
    real(dp), intent(in) :: nu
    real(dp)             :: f2t, y

    f2t = 1.d0/(compute_bnu_prime(nu))*1.0d-14

  end function f2t

  subroutine convert_maps(self,dat)
    ! We want to run everything in uK_RJ, yeah?
    implicit none
    type(params), intent(inout) :: self
    type(data),   intent(inout) :: dat
    integer(i4b)                :: j
    
    do j = 1, nbands
       if (.not. self%bp_map(j)) then
          if (trim(self%dat_unit(j)) == 'uK_RJ') then
             cycle
          else if (trim(self%dat_unit(j)) == 'uK_cmb') then
             ! uK_cmb -> uK_RJ
             write(*,*) 'Putting band ', trim(self%dat_label(j)), ' from uK_cmb to uK_RJ.'
             dat%sig_map(:,:,j) = dat%sig_map(:,:,j)/a2t(self%dat_nu(j)*1.0d9)
             dat%rms_map(:,:,j) = dat%rms_map(:,:,j)/a2t(self%dat_nu(j)*1.0d9)
          else if (trim(self%dat_unit(j)) == 'MJy/sr') then
             ! MJy/sr -> uK_RJ
             write(*,*) 'Putting band ', trim(self%dat_label(j)), ' from MJy/sr to uK_RJ'
             dat%sig_map(:,:,j) = dat%sig_map(:,:,j)/a2f(self%dat_nu(j)*1.0d9)
             dat%rms_map(:,:,j) = dat%rms_map(:,:,j)/a2f(self%dat_nu(j)*1.0d9)
          else
             write(*,*) 'Not a unit, dumbass!'
             stop
          end if
       end if
    end do
    
  end subroutine convert_maps
  
  subroutine convert_maps_bp(self,dat)
    implicit none
    type(params), intent(inout) :: self
    type(data),   intent(inout) :: dat
    integer(i4b)                :: j
    
    do j = 1, nbands
       if (self%bp_map(j)) then
          if (trim(self%dat_unit(j)) == 'uK_RJ') then
             cycle
          else if (trim(self%dat_unit(j)) == 'uK_cmb') then
             ! uK_cmb -> uK_RJ
             write(*,*) 'Putting band ', trim(self%dat_label(j)), ' from uK_cmb to uK_RJ.'
             dat%sig_map(:,:,j) = dat%sig_map(:,:,j)/a2t(self%dat_nu(j)*1.0d9)
             dat%rms_map(:,:,j) = dat%rms_map(:,:,j)/a2t(self%dat_nu(j)*1.0d9)
          else if (trim(self%dat_unit(j)) == 'MJy/sr') then
             ! MJy/sr -> uK_RJ
             write(*,*) 'Putting band ', trim(self%dat_label(j)), ' from MJy/sr to uK_RJ'
             dat%sig_map(:,:,j) = dat%sig_map(:,:,j)/a2f(self%dat_nu(j)*1.0d9)
             dat%rms_map(:,:,j) = dat%rms_map(:,:,j)/a2f(self%dat_nu(j)*1.0d9)
          else
             write(*,*) 'Not a unit, dumbass!'
             stop
          end if
       end if
    end do
    
  end subroutine convert_maps_bp


end module dang_data_mod
