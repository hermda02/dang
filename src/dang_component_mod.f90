module dang_component_mod
  use healpix_types
  use pix_tools
  use fitstools
  use dang_util_mod
  use dang_param_mod
  use dang_bp_mod
  implicit none
  
  public dang_comps, component_pointer, component_list

  type :: dang_comps
     
     character(len=16)                                :: label, type
     logical(lgt)                                     :: sample_amplitude
     logical(lgt),      allocatable, dimension(:)     :: sample_index
     real(dp)                                         :: nu_ref
     integer(i4b)                                     :: cg_group
     logical(lgt)                                     :: polfit
     logical(lgt),      allocatable, dimension(:)     :: corr             ! Do we correct this band?
     integer(i4b)                                     :: nfit

     real(dp),          allocatable, dimension(:,:,:) :: indices
     real(dp),          allocatable, dimension(:,:)   :: amplitude

     real(dp),          allocatable, dimension(:,:)   :: spec_par
     real(dp),          allocatable, dimension(:)     :: temp_amps

     real(dp),          allocatable, dimension(:,:)   :: beta_s, beta_d, T_d, HI
     real(dp),          allocatable, dimension(:)     :: HI_amps

     character(len=16), allocatable, dimension(:)     :: prior_type
     real(dp),          allocatable, dimension(:,:)   :: gauss_prior
     real(dp),          allocatable, dimension(:,:)   :: uni_prior
    
   contains

     procedure :: eval_sed

  end type dang_comps

  interface dang_comps
     procedure constructor
  end interface dang_comps

  type component_pointer
     type(dang_comps), pointer :: p => null()
  end type component_pointer

  type(component_pointer), allocatable, dimension(:) :: component_list 

contains 

  function constructor(dpar,component)
    ! The overall constructor for the component class
    ! Takes the parameters as read in from the parameter 
    ! file and sets them as attributes to the class
    !
    ! Input: 
    !        dpar: class - dang_params
    !        component: integer - which component are we initializing?
    !    
    implicit none

    type(dang_params), intent(in) :: dpar
    class(dang_comps), pointer    :: constructor
    integer(i4b),      intent(in) :: component

    allocate(constructor)

    constructor%label            = trim(dpar%fg_label(component))
    constructor%type             = dpar%fg_type(component)
    constructor%sample_amplitude = dpar%fg_samp_amp(component)
    constructor%cg_group         = dpar%fg_cg_group(component)

    if (trim(constructor%type) == 'mbb') then
       ! Allocate arrays to appropriate size for each component type
       allocate(constructor%gauss_prior(2,2))
       allocate(constructor%uni_prior(2,2))
       allocate(constructor%sample_index(2))
       allocate(constructor%prior_type(2))

       ! Allocate maps for the components
       allocate(constructor%amplitude(0:npix-1,nmaps))
       allocate(constructor%indices(0:npix-1,nmaps,2))
     
       ! Bools for sampling indices
       constructor%sample_index(1)  = dpar%fg_samp_spec(component,1)
       constructor%sample_index(2)  = dpar%fg_samp_spec(component,2)

       ! Define prior for likelihood evaluation
       constructor%prior_type(1)    = dpar%fg_prior_type(component,1)
       constructor%prior_type(2)    = dpar%fg_prior_type(component,2)

       ! Define the Gaussian prior and bounds
       ! Beta
       constructor%gauss_prior(1,1) = dpar%fg_gauss(component,1,1)
       constructor%gauss_prior(1,2) = dpar%fg_gauss(component,1,2)
       ! T_d
       constructor%gauss_prior(2,1) = dpar%fg_gauss(component,2,1)
       constructor%gauss_prior(2,2) = dpar%fg_gauss(component,2,2)
       
       ! Beta
       constructor%uni_prior(1,1)   = dpar%fg_uni(component,1,1)
       constructor%uni_prior(1,2)   = dpar%fg_uni(component,1,2)
       ! T_d
       constructor%uni_prior(2,1)   = dpar%fg_uni(component,2,1)
       constructor%uni_prior(2,2)   = dpar%fg_uni(component,2,2)

       constructor%nu_ref           = dpar%fg_nu_ref(component) 

    else if (trim(constructor%type) == 'power-law') then
       ! Allocate arrays to appropriate size for each component type
       allocate(constructor%gauss_prior(1,2))
       allocate(constructor%uni_prior(1,2))
       allocate(constructor%sample_index(1))
       allocate(constructor%prior_type(1))

       ! Allocate maps for the components
       allocate(constructor%amplitude(0:npix-1,nmaps))
       allocate(constructor%indices(0:npix-1,nmaps,1))
       
       ! Bools for sampling indices
       constructor%sample_index(1)  = dpar%fg_samp_spec(component,1)

       ! Define prior for likelihood evaluation
       constructor%prior_type(1)    = dpar%fg_prior_type(component,1)
       
       ! Define the Gaussian prior and bounds
       constructor%gauss_prior(1,1) = dpar%fg_gauss(component,1,1)
       constructor%gauss_prior(1,2) = dpar%fg_gauss(component,1,2)

       constructor%uni_prior(1,1)   = dpar%fg_uni(component,1,1)
       constructor%uni_prior(1,2)   = dpar%fg_uni(component,1,2)

       constructor%nu_ref           = dpar%fg_nu_ref(component) 

       constructor%indices          = dpar%fg_init(1,1)

    else if (trim(constructor%type) == 'template') then
       
       allocate(constructor%corr(nbands))

       constructor%polfit           = .true. ! WARNING HARD CODED TO .true. FOR NOW
       constructor%nfit             = dpar%fg_nfit(component) ! Also currently hardcoded
       constructor%corr             = dpar%fg_temp_corr(component,:)

    else
       write(*,*) "Warning - unrecognized component type detected"
       stop

    end if
  end function constructor

  subroutine initialize_components(dpar)
    implicit none
    type(dang_params)             :: dpar
    integer(i4b)                  :: i

    allocate(component_list(dpar%ncomp))

    do i = 1, dpar%ncomp
       write(*,*) 'Initialize component ', i
       component_list(i)%p => dang_comps(dpar,i)
    end do

  end subroutine initialize_components

  subroutine init_synch(self,dpar,npix,nmaps)
    implicit none
    type(dang_comps)              :: self
    type(dang_params)             :: dpar
    integer(i4b), intent(in)      :: npix, nmaps
    
    allocate(self%beta_s(0:npix-1,nmaps))
    write(*,*) 'Allocated synch parameter maps'
    if (trim(dpar%fg_spec_file(1,1)) == 'none') then 
       write(*,fmt='(a,f8.4)') 'Full sky beta_s estimate ', dpar%fg_init(1,1)
       write(*,*) ''
       self%beta_s     = dpar%fg_init(1,1) ! Synchrotron beta initial guess
    else
       write(*,*) "No init files found"
       stop
       !call read_bintab(trim(param%fg_spec_file(1,1)),self%beta_s,npix,3,nullval,anynull,header=header)
    end if
    
  end subroutine init_synch
  
  subroutine init_dust(self,dpar,npix,nmaps)
    implicit none
    type(dang_comps)              :: self
    type(dang_params)             :: dpar
    integer(i4b), intent(in)      :: npix, nmaps
    
    allocate(self%beta_d(0:npix-1,nmaps))
    allocate(self%T_d(0:npix-1,nmaps))
    write(*,*) 'Allocated dust parameter maps!'
    write(*,*) ''
    self%beta_d     = 1.53d0              ! Dust beta initial guess
    self%T_d        = 19.6d0
    
  end subroutine init_dust

  subroutine init_hi_fit(self, dpar, npix)
    implicit none
    type(dang_comps)         :: self
    type(dang_params)        :: dpar
    integer(i4b), intent(in) :: npix
    character(len=80), dimension(180) :: head

    allocate(self%HI(0:npix-1,1))
    allocate(self%T_d(0:npix-1,nmaps))
    allocate(self%HI_amps(dpar%numinc))
    write(*,*) 'Allocated HI fitting maps!'
    write(*,*) ''
    call read_bintab(trim(dpar%datadir)//trim(dpar%HI_file),self%HI,npix,1,nullval,anynull,header=head)

    if (trim(dpar%HI_Td_init) == 'none') then
       self%T_d = dpar%HI_Td_mean
    else
       call read_bintab(trim(dpar%datadir)//trim(dpar%HI_Td_init),self%T_d,npix,1,nullval,anynull,header=head)
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
  ! function compute_spectrum(dpar, self, bp, ind, freq, pix, map_n, index)

  function eval_sed(self, dpar, bp, ind, pix, map_n, index)
    ! always computed in RJ units
    
    implicit none
    class(dang_params)                 :: dpar
    class(dang_comps)                  :: self
    type(bandinfo)                     :: bp
    integer(i4b),           intent(in) :: ind
    integer(i4b),           optional   :: pix
    integer(i4b),           optional   :: map_n
    integer(i4b)                       :: i
    real(dp), dimension(:), optional   :: index
    real(dp)                           :: z, compute_spectrum
    real(dp)                           :: eval_sed

    if (trim(self%type) == 'power-law') then
       if (bp%id == 'delta') then
          ! if (ind == 1) then
          if (present(index)) then
             compute_spectrum = (bp%nu_c/self%nu_ref)**index(1)
          else 
             compute_spectrum = (bp%nu_c/self%nu_ref)**self%beta_s(pix,map_n)
          end if
       else
          compute_spectrum = 0.d0
          ! Compute for LFI bandpass
          if (present(index)) then
             do i = 1, bp%n
                compute_spectrum = compute_spectrum + bp%tau0(i)*(bp%nu0(i)/self%nu_ref)**index(1)
             end do
          else 
             do i = 1, bp%n
                compute_spectrum = compute_spectrum + bp%tau0(i)*(bp%nu0(i)/self%nu_ref)**self%indices(pix,map_n,1)
             end do
          end if
       end if
    else if (trim(self%type) == 'mbb') then
       if (bp%id == 'delta') then
          if (present(index)) then
             z = h / (k_B*index(2))
             compute_spectrum = (exp(z*self%nu_ref*1d9)-1.d0) / &
                  & (exp(z*bp%nu_c)-1.d0) * (bp%nu_c/(self%nu_ref*1d9))**(index(1)+1.d0)
          else
             z = h / (k_B*self%T_d(pix,map_n))
             compute_spectrum = (exp(z*self%nu_ref*1d9)-1.d0) / &
                  & (exp(z*bp%nu_c)-1.d0) * (bp%nu_c/(self%nu_ref*1d9))**(self%beta_d(pix,map_n)+1.d0)
          end if
       else
          compute_spectrum = 0.d0
          if (present(index)) then
             z = h / (k_B*index(2))
             do i = 1, bp%n
                compute_spectrum = compute_spectrum + bp%tau0(i)*(exp(z*self%nu_ref*1d9)-1.d0) / &
                     (exp(z*bp%nu0(i))-1.d0) * (bp%nu0(i)/(self%nu_ref*1d9))**(index(1)+1.d0)
             end do
          else
             z = h / (k_B*self%T_d(pix,map_n))
             do i = 1, bp%n
                compute_spectrum = compute_spectrum + bp%tau0(i)*(exp(z*self%nu_ref*1d9)-1.d0) / &
                     (exp(z*bp%nu0(i))-1.d0) * (bp%nu0(i)/(self%nu_ref*1d9))**(self%beta_d(pix,map_n)+1.d0)
             end do
          end if
       end if
    end if
    eval_sed = compute_spectrum

  end function eval_sed
  
  ! function compute_spectrum(dpar, self, bp, ind, freq, pix, map_n, index)
  function compute_spectrum(dpar, self, bp, ind, pix, map_n, index)
    ! always computed in RJ units
    
    implicit none
    class(dang_params)             :: dpar
    type(dang_comps)               :: self
    type(bandinfo)                 :: bp
    integer(i4b),       intent(in) :: ind
    integer(i4b),       optional   :: pix
    integer(i4b),       optional   :: map_n
    integer(i4b)                   :: i
    real(dp),           optional   :: index
    real(dp)                       :: z, compute_spectrum

    ! if (trim(dpar%fg_label(ind)) == 'power-law') then
    if (bp%id == 'delta') then
       if (ind == 1) then
          if (present(index)) then
             compute_spectrum = (bp%nu_c/dpar%fg_nu_ref(ind))**index
          else 
             compute_spectrum = (bp%nu_c/dpar%fg_nu_ref(ind))**self%beta_s(pix,map_n)
          end if
          !else if (trim(dpar%fg_label(ind)) == 'mbb') then
       else if (ind == 2) then
          write(*,*) 'yuh'
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
                compute_spectrum = compute_spectrum + bp%tau0(i)*(bp%nu0(i)/dpar%fg_nu_ref(ind))**index
             end do
          else 
             do i = 1, bp%n
                compute_spectrum = compute_spectrum + bp%tau0(i)*(bp%nu0(i)/dpar%fg_nu_ref(ind))**self%beta_s(pix,map_n)
             end do
          end if
       ! else if (trim(dpar%fg_label(ind)) == 'mbb') then
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
