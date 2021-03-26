module dang_component_mod
  use healpix_types
  use pix_tools
  use fitstools
  use dang_util_mod
  use dang_param_mod
  implicit none
  
  type, public                              :: component
     
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
       write(*,*) 'Full sky beta_s estimate ', param%fg_init(1,1)
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
  
end module dang_component_mod
