module dang_component_mod
  use healpix_types
  use pix_tools
  use fitstools
  use dang_util_mod
  use dang_param_mod
  use dang_data_mod
  implicit none
  
  type, public                              :: component
     
     character(len=32), allocatable, dimension(:) :: joint
     real(dp), allocatable, dimension(:,:)        :: beta_s, beta_d, T_d, HI
     real(dp), allocatable, dimension(:)          :: HI_amps

  end type component

contains 

  subroutine init_synch(self,param,npix,nmaps)
    implicit none
    type(component)               :: self
    type(params)                  :: param
    integer(i4b), intent(in)      :: npix, nmaps
    
    allocate(self%beta_s(0:npix-1,nmaps))
    write(*,*) 'Allocated synch maps'
    if (trim(param%fg_spec_file(1,1)) == 'none') then 
       self%beta_s     = param%fg_init(1,1) ! Synchrotron beta initial guess
    else
       !call read_bintab(trim(param%fg_spec_file(1,1)),self%beta_s,npix,3,nullval,anynull,header=header)
    end if
    
  end subroutine init_synch
  
  subroutine init_dust(self,param,npix,nmaps)
    implicit none
    type(component)               :: self
    type(params)                  :: param
    integer(i4b), intent(in)      :: npix, nmaps
    
    allocate(self%beta_d(0:npix-1,nmaps))
    allocate(self%T_d(0:npix-1,nmaps))
    write(*,*) 'Allocated dust maps!'

    self%beta_d     = 1.53d0              ! Dust beta initial guess
    self%T_d        = 19.6d0
    
  end subroutine init_dust

  subroutine init_template(self,npix,nmaps,ntemp,nbands)
    implicit none
    type(data)               :: self
    integer(i4b), intent(in) :: npix, nmaps, ntemp, nbands

    allocate(self%temps(0:npix-1,nmaps,ntemp))
    allocate(self%temp_amps(nbands,nmaps,ntemp))
    allocate(self%temp_norm(nmaps,ntemp))

  end subroutine init_template

  subroutine init_hi_fit(self, param, npix)
    implicit none
    type(component)          :: self
    type(params)             :: param
    integer(i4b), intent(in) :: npix
    character(len=80), dimension(180) :: head

    allocate(self%HI(0:npix-1,1))
    allocate(self%T_d(0:npix-1,nmaps))
    allocate(self%HI_amps(param%numinc))
    write(*,*) 'Allocated HI fitting maps!'

    call read_bintab(trim(param%datadir)//trim(param%HI_file),self%HI,npix,1,nullval,anynull,header=head)

    if (trim(param%HI_Td_init) == 'none') then
       self%T_d = param%HI_Td_mean
    else
       call read_bintab(trim(param%HI_Td_init),self%T_d,npix,1,nullval,anynull,header=head)
    end if

  end subroutine init_hi_fit
  
  function planck(fre,T)
    implicit none
    real(dp), intent(in)  :: fre
    real(dp), intent(in)  :: T
    real(dp)              :: planck
    ! Output in units of [W sr^-1 m^-2 Hz^-1]
    planck  = ((2.d0*h*fre**3.d0)/(c**2.d0))*(1.d0/(exp((h*fre)/(k_B*T))-1))
  end function planck
  
  function compute_spectrum(param, self, ind, freq, pix, mapn, index)
    ! always computed in RJ units
    
    implicit none
    class(params)                  :: param
    type(component)                :: self
    real(dp),           intent(in) :: freq
    integer(i4b),       intent(in) :: ind
    integer(i4b),       intent(in) :: pix
    integer(i4b),       intent(in) :: mapn
    real(dp), optional             :: index
    real(dp)                       :: z, compute_spectrum
    
    !if (trim(param%fg_label(ind)) == 'power-law') then
    if (ind == 1) then
      if (present(index)) then
          compute_spectrum = (freq/param%fg_nu_ref(ind))**index
       else 
          compute_spectrum = (freq/param%fg_nu_ref(ind))**self%beta_s(pix,mapn)
       end if
    !else if (trim(param%fg_label(ind)) == 'mbb') then
    else if (ind == 2) then
       z = h / (k_B*self%T_d(pix,mapn))
!       compute_spectrum = (exp(z*param%fg_nu_ref(ind)*1d9)-1.d0) / &
!            (exp(z*freq*1d9)-1.d0) * (freq/param%fg_nu_ref(ind))**(self%beta_d(pix,mapn)+1.d0)!*rj_cmb
       compute_spectrum = (exp(z*353.d0*1d9)-1.d0) / &
            (exp(z*freq*1d9)-1.d0) * (freq/353.d0)**(self%beta_d(pix,mapn)+1.d0)
    end if
  end function compute_spectrum

  subroutine dust_correct_band(dat,param,comp,band)
    implicit none
    type(data),   intent(inout) :: dat
    type(params)                :: param
    type(component)             :: comp
    integer(i4b), intent(in)    :: band
    real(dp), allocatable, dimension(:,:,:) :: thermal_map
    integer(i4b)                :: i, j, k
    character(len=256)          :: title

    allocate(thermal_map(0:npix-1,nmaps,nbands))

    if (trim(param%dust_corr_type) == 'uniform') then
       comp%T_d    = param%mbb_gauss(1,1)
       comp%beta_d = param%mbb_gauss(2,1)
    else if (trim(param%dust_corr_type) == 'sample') then
       if (param%mbb_gauss(1,2) .gt. 0.d0) then
          comp%T_d    = rand_normal(param%mbb_gauss(1,1),param%mbb_gauss(1,2))
       else 
          comp%T_d    = param%mbb_gauss(1,1)
       end if
       if (param%mbb_gauss(2,2) .gt. 0.d0) then
          comp%beta_d = rand_normal(param%mbb_gauss(2,1),param%mbb_gauss(2,2))
       else
          comp%beta_d = param%mbb_gauss(2,1)
       end if
    else if (trim(param%dust_corr_type) == 'planck') then
       stop
    end if
    write(*,'(a,a)') 'Dust correcting band ', trim(param%dat_label(band))
    do k = param%pol_type(1), param%pol_type(size(param%pol_type))
       do i = 0, npix-1
          thermal_map(i,k,band) = dat%temps(i,k,1)*compute_spectrum(param,comp,2,param%dat_nu(band),i,k)
          dat%sig_map(i,k,band) = dat%sig_map(i,k,band) - thermal_map(i,k,band)
       end do
    end do
    title = trim(param%outdir)//trim(param%dat_label(band))//'_thermal_map.fits'
    call write_result_map(trim(title), nside, ordering, header, thermal_map(:,:,band))
  end subroutine dust_correct_band

  subroutine extrapolate_foreground(param, dat, comp, ind, map_n)
    implicit none
    type(data),   intent(inout) :: dat
    type(params)                :: param
    type(component)             :: comp
    integer(i4b), intent(in)    :: ind, map_n
    integer(i4b)                :: i, j, k

    do i = 0, npix-1
       do j = 1, nbands
          do k = param%pol_type(1), param%pol_type(size(param%pol_type))
             dat%fg_map(i,k,j,ind) = dat%fg_map(i,k,param%fg_ref_loc(ind),ind)*compute_spectrum(param,comp,ind,param%dat_nu(j),i,k)
          end do
       end do
    end do

  end subroutine extrapolate_foreground

  subroutine extrapolate_template(param, dat, comp, ind, map_n)
    implicit none
    type(data),   intent(inout) :: dat
    type(params)                :: param
    type(component)             :: comp
    integer(i4b), intent(in)    :: ind, map_n
    integer(i4b)                :: i, j, k

    do i = 0, npix-1
       do j = 1, nbands
          do k = param%pol_type(1), param%pol_type(size(param%pol_type))
             dat%fg_map(i,k,j,param%ncomp+ind) = dat%temp_amps(j,k,ind)*dat%temps(i,k,ind)
          end do
       end do
    end do

  end subroutine extrapolate_template
  
end module dang_component_mod
