module dang_component_mod
  use healpix_types
  use utility_mod
  use dang_param_mod
  use dang_data_mod
  implicit none
  
  type, public                              :: component
     
     character(len=32), allocatable, dimension(:) :: joint
     real(dp), allocatable, dimension(:,:)        :: beta_s, beta_d, T_d

  end type component

contains 

  subroutine init_synch(self,npix,nmaps)
    implicit none
    type(component)               :: self
    integer(i4b), intent(in) :: npix, nmaps
    
    allocate(self%beta_s(0:npix-1,nmaps))
    
  end subroutine init_synch
  
  subroutine init_dust(self,npix,nmaps)
    implicit none
    type(component)               :: self
    integer(i4b), intent(in) :: npix, nmaps
    
    allocate(self%beta_d(0:npix-1,nmaps))
    allocate(self%T_d(0:npix-1,nmaps))
    
  end subroutine init_dust
  
  function planck(fre,T)
    implicit none
    real(dp), intent(in)  :: fre
    real(dp), intent(in)  :: T
    real(dp)              :: planck
    ! Output in units of [W sr^-1 m^-2 Hz^-1]
    planck  = ((2.d0*h*fre**3.d0)/(c**2.d0))*(1.d0/(exp((h*fre)/(k_B*T))-1))
  end function planck
  
  function compute_spectrum(self, comp, ind, freq, pix, mapn, param)
    ! always computed in RJ units
    
    implicit none
    class(params)                  :: self
    type(component)                :: comp
    real(dp),           intent(in) :: freq
    integer(i4b),       intent(in) :: ind
    integer(i4b),       intent(in) :: pix
    integer(i4b),       intent(in) :: mapn
    real(dp), optional             :: param
    real(dp)                       :: z, compute_spectrum
    
    if (ind == 1) then
       if (present(param)) then
          compute_spectrum = (freq/self%fg_nu_ref(ind))**param
       else 
          compute_spectrum = (freq/self%fg_nu_ref(ind))**comp%beta_s(pix,mapn)
       end if
!    else if (trim(self%fg_label(ind)) == 'mbb') then
    else if (ind == 2) then
       z = h / (k_B*comp%T_d(pix,mapn))
!       compute_spectrum = (exp(z*self%fg_nu_ref(ind)*1d9)-1.d0) / &
!            (exp(z*freq*1d9)-1.d0) * (freq/self%fg_nu_ref(ind))**(comp%beta_d(pix,mapn)+1.d0)!*rj_cmb
       compute_spectrum = (exp(z*353.d0*1d9)-1.d0) / &
            (exp(z*freq*1d9)-1.d0) * (freq/353.d0)**(comp%beta_d(pix,mapn)+1.d0)
    end if
  end function compute_spectrum

  subroutine dust_correct_band(dat,param,comp,band)
    implicit none
    type(data),   intent(inout) :: dat
    type(params)                :: param
    type(component)             :: comp
    integer(i4b), intent(in)    :: band

    integer(i4b)                :: i, j, k

    if (trim(param%dust_corr_type) == 'uniform') then
       comp%beta_d = 1.53d0
       comp%T_d    = 19.6d0
    else if (trim(param%dust_corr_type) == 'sample') then
       comp%beta_d = random_number(1.53d0,0.02d0)
       comp%T_d    = 19.6d0
    else if (trim(param%dust_corr_type) == 'planck') then
       stop
       ! call read_bintab()
       ! call read_bintab()
    end if
    write(*,'(a,a)') 'Dust correcting band ', trim(param%dat_label(band))
    do k = param%pol_type(1), param%pol_type(size(param%pol_type))
       do i = 0, npix-1
          dat%sig_map(i,k,band) = dat%sig_map(i,k,band) - dat%temps(i,k,1)*compute_spectrum(param,comp,2,param%dat_nu(band),i,k)
       end do
    end do
  end subroutine dust_correct_band
  
end module dang_component_mod
