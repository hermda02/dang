module dang_component_mod
  use healpix_types
  use pix_tools
  use fitstools
  use dang_util_mod
  use dang_param_mod
  use dang_bp_mod
  implicit none
  
  type, public                              :: dang_comps
     
     real(dp), allocatable, dimension(:,:)        :: beta_s, beta_d, T_d, HI
     real(dp), allocatable, dimension(:)          :: HI_amps

  end type dang_comps

contains 

  subroutine init_synch(self,param,npix,nmaps)
    implicit none
    type(dang_comps)              :: self
    type(dang_params)             :: param
    integer(i4b), intent(in)      :: npix, nmaps
    
    allocate(self%beta_s(0:npix-1,nmaps))
    write(*,*) 'Allocated synch parameter maps'
    if (trim(param%fg_spec_file(1,1)) == 'none') then 
       write(*,fmt='(a,f8.4)') 'Full sky beta_s estimate ', param%fg_init(1,1)
       self%beta_s     = param%fg_init(1,1) ! Synchrotron beta initial guess
    else
       write(*,*) "No init files found"
       stop
       !call read_bintab(trim(param%fg_spec_file(1,1)),self%beta_s,npix,3,nullval,anynull,header=header)
    end if
    
  end subroutine init_synch
  
  subroutine init_dust(self,param,npix,nmaps)
    implicit none
    type(dang_comps)              :: self
    type(dang_params)             :: param
    integer(i4b), intent(in)      :: npix, nmaps
    
    allocate(self%beta_d(0:npix-1,nmaps))
    allocate(self%T_d(0:npix-1,nmaps))
    write(*,*) 'Allocated dust parameter maps!'

    self%beta_d     = 1.53d0              ! Dust beta initial guess
    self%T_d        = 19.6d0
    
  end subroutine init_dust

  subroutine init_hi_fit(self, param, npix)
    implicit none
    type(dang_comps)         :: self
    type(dang_params)        :: param
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
  
  function planck(bp,T)
    implicit none
    type(bandinfo), intent(in) :: bp
    real(dp),       intent(in) :: T
    real(dp)                   :: planck
    integer(i4b)               :: i
    ! Output in units of [W sr^-1 m^-2 Hz^-1]
    if (bp%id == 'delta') then
       planck  = ((2.d0*h*(bp%nu_c)**3.d0)/(c**2.d0))*(1.d0/(exp((h*bp%nu_c)/(k_B*T))-1))
    else
       planck = 0.d0
       do i = 1, bp%n
          planck = planck + bp%tau0(i)*((2.d0*h*(bp%nu0(i))**3.d0)/(c**2.d0))*(1.d0/(exp((h*(bp%nu0(i)))/(k_B*T))-1))
       end do
    end if
  end function planck
  
  ! function compute_spectrum(param, self, bp, ind, freq, pix, map_n, index)
  function compute_spectrum(param, self, bp, ind, pix, map_n, index)
    ! always computed in RJ units
    
    implicit none
    class(dang_params)             :: param
    type(dang_comps)               :: self
    type(bandinfo)                 :: bp
    integer(i4b),       intent(in) :: ind
    integer(i4b),       optional   :: pix
    integer(i4b),       optional   :: map_n
    integer(i4b)                   :: i
    real(dp),           optional   :: index
    real(dp)                       :: z, compute_spectrum
    
    !if (trim(param%fg_label(ind)) == 'power-law') then
    if (bp%id == 'delta') then
       if (ind == 1) then
          if (present(index)) then
             compute_spectrum = (bp%nu_c/param%fg_nu_ref(ind))**index
          else 
             compute_spectrum = (bp%nu_c/param%fg_nu_ref(ind))**self%beta_s(pix,map_n)
          end if
          !else if (trim(param%fg_label(ind)) == 'mbb') then
       else if (ind == 2) then
          z = h / (k_B*self%T_d(pix,map_n))
          compute_spectrum = (exp(z*353.d0*1d9)-1.d0) / &
               (exp(z*bp%nu_c)-1.d0) * (bp%nu_c/(353.d0*1d9))**(self%beta_d(pix,map_n)+1.d0)
       end if
    else
       compute_spectrum = 0.d0
       ! Compute for LFI bandpass
       if (ind == 1) then
          if (present(index)) then
             do i = 1, bp%n
                compute_spectrum = compute_spectrum + bp%tau0(i)*(bp%nu0(i)/param%fg_nu_ref(ind))**index
             end do
          else 
             do i = 1, bp%n
                compute_spectrum = compute_spectrum + bp%tau0(i)*(bp%nu0(i)/param%fg_nu_ref(ind))**self%beta_s(pix,map_n)
             end do
          end if
          !else if (trim(param%fg_label(ind)) == 'mbb') then
       else if (ind == 2) then
          z = h / (k_B*self%T_d(pix,map_n))
          do i = 1, bp%n
             compute_spectrum = compute_spectrum + bp%tau0(i)*(exp(z*353.d0*1d9)-1.d0) / &
               (exp(z*bp%nu0(i))-1.d0) * (bp%nu0(i)/353.d9)**(self%beta_d(pix,map_n)+1.d0)
          end do
       end if
    end if
  end function compute_spectrum
  
end module dang_component_mod
