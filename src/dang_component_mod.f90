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

     ! General component information
     integer(i4b)                                     :: cg_group         ! CG group number 
     character(len=16)                                :: label, type      ! Component label and type
     integer(i4b)                                     :: nfit             ! (Template only)
     integer(i4b)                                     :: nindices         ! How many indices does the component have?
     real(dp)                                         :: nu_ref           ! Reference frequency
     logical(lgt)                                     :: sample_amplitude ! Do we sample this spectral index?

     logical(lgt),      allocatable, dimension(:)     :: corr             ! Do we correct this band?
     character(len=16), allocatable, dimension(:)     :: ind_label        ! Label of the indices
     logical(lgt),      allocatable, dimension(:)     :: sample_index     ! Do we sample this spectral index?
     real(dp),          allocatable, dimension(:)     :: step_size        ! M-H step size

     ! Sky/template information
     real(dp),          allocatable, dimension(:,:)   :: amplitude           ! Amplitude maps
     real(dp),          allocatable, dimension(:,:)   :: template            ! Template map
     real(dp),          allocatable, dimension(:,:)   :: template_amplitudes ! Template amplitudes

     integer(i4b),      allocatable, dimension(:)     :: index_mode          ! Fullsky/per-pixel
     real(dp),          allocatable, dimension(:,:,:) :: indices             ! Indices maps

     ! Prior information
     character(len=16), allocatable, dimension(:)     :: lnl_type
     character(len=16), allocatable, dimension(:)     :: prior_type
     real(dp),          allocatable, dimension(:,:)   :: gauss_prior
     real(dp),          allocatable, dimension(:,:)   :: uni_prior

     ! Poltype information
     integer(i4b),      allocatable, dimension(:)     :: nflag
     integer(i4b),      allocatable, dimension(:,:)   :: pol_flag
    
   contains

     procedure :: eval_sed
     procedure :: eval_signal

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

    integer(i4b)                  :: i, j, count

    integer(i4b), allocatable, dimension(:) :: flag_buffer

    allocate(constructor)

    constructor%label            = trim(dpar%fg_label(component))
    constructor%type             = dpar%fg_type(component)
    constructor%cg_group         = dpar%fg_cg_group(component)

    if (trim(constructor%type) /= 'template' .and. trim(constructor%type) /= 'hi_fit') then
       constructor%sample_amplitude = dpar%fg_amp_samp(component)
    else
       constructor%sample_amplitude = .true.
    end if


    if (trim(constructor%type) == 'mbb') then

       ! Allocate arrays to appropriate size for each component type
       constructor%nindices = 2
       allocate(constructor%gauss_prior(constructor%nindices,2))
       allocate(constructor%uni_prior(constructor%nindices,2))
       allocate(constructor%sample_index(constructor%nindices))
       allocate(constructor%prior_type(constructor%nindices))
       allocate(constructor%lnl_type(constructor%nindices))
       allocate(constructor%index_mode(constructor%nindices))
       allocate(constructor%step_size(constructor%nindices))

       ! Allocate general pol_type flag array
       allocate(constructor%nflag(constructor%nindices))
       allocate(constructor%pol_flag(constructor%nindices,3)) ! The three is here because that's the max number of poltypes we can handle

       ! Allocate maps for the components
       allocate(constructor%amplitude(0:npix-1,nmaps))
       allocate(constructor%indices(0:npix-1,nmaps,constructor%nindices))

       allocate(constructor%ind_label(constructor%nindices))

       constructor%ind_label = ['BETA', 'T']

       ! Reference frequency
       constructor%nu_ref           = dpar%fg_nu_ref(component) 

       ! Initialize an amplitude map, or don't
       if (trim(dpar%fg_filename(component)) == 'none') then
          constructor%amplitude = 0.d0
       else
          call read_bintab(trim(dpar%datadir)//trim(dpar%fg_filename(component)),&
               constructor%amplitude, npix, nmaps, nullval, anynull, header=header)
       end if

       ! Initialize pol_flag arrays
       constructor%nflag = 0
       constructor%pol_flag = 0

       do i = 1, constructor%nindices

          ! Load the poltype into the flag buffer, and store bit flags for each component
          flag_buffer = return_poltype_flag(dpar%fg_spec_poltype(component,i))
          constructor%nflag = size(flag_buffer)
          do j = 1, size(flag_buffer)
             constructor%pol_flag(i,j) = flag_buffer(j)
          end do

          ! Do we sample this index?
          constructor%sample_index(i)  = dpar%fg_samp_spec(component,i)

          ! Sample full sky or per-pixel?
          if (trim(dpar%fg_ind_region(component,i)) == 'fullsky') then
             constructor%index_mode(i) = 1
          else if (trim(dpar%fg_ind_region(component,i)) == 'per-pixel') then
             constructor%index_mode(i) = 2
          end if

          ! Define the lnl evaluation for each index
          constructor%lnl_type(i)      = dpar%fg_ind_lnl(component,i)

          ! Define MH step siz
          constructor%step_size(i) = dpar%fg_spec_step(component,i)

          ! Define prior for likelihood evaluation
          constructor%prior_type(i)    = dpar%fg_prior_type(component,i) 
          constructor%gauss_prior(i,1) = dpar%fg_gauss(component,i,1)
          constructor%gauss_prior(i,2) = dpar%fg_gauss(component,i,2)
          constructor%uni_prior(i,1)   = dpar%fg_uni(component,i,1)
          constructor%uni_prior(i,2)   = dpar%fg_uni(component,i,2)

          ! Initialize spectral index maps, or don't
          if (trim(dpar%fg_spec_file(component,i)) == 'none') then
             constructor%indices(:,:,i)   = dpar%fg_init(component,i)
          else
             call read_bintab(trim(dpar%datadir)//trim(dpar%fg_spec_file(component,i)),&
                  constructor%indices(:,:,i), npix, nmaps, nullval, anynull, header=header)
          end if
       end do

    else if (trim(constructor%type) == 'power-law') then
       ! Allocate arrays to appropriate size for each component type
       constructor%nindices = 1
       allocate(constructor%gauss_prior(constructor%nindices,2))
       allocate(constructor%uni_prior(constructor%nindices,2))
       allocate(constructor%sample_index(constructor%nindices))
       allocate(constructor%prior_type(constructor%nindices))
       allocate(constructor%lnl_type(constructor%nindices))
       allocate(constructor%index_mode(constructor%nindices))
       allocate(constructor%step_size(constructor%nindices))

       ! Allocate general pol_type flag array
       allocate(constructor%nflag(constructor%nindices))
       allocate(constructor%pol_flag(constructor%nindices,3)) ! The three is here because that's the max number of poltypes we can handle
       
       ! Allocate maps for the components
       allocate(constructor%amplitude(0:npix-1,nmaps))
       allocate(constructor%indices(0:npix-1,nmaps,constructor%nindices))
       allocate(constructor%ind_label(constructor%nindices))

       constructor%ind_label = ['T']
       
       ! Reference frequency
       constructor%nu_ref           = dpar%fg_nu_ref(component) 

       ! Initialize an amplitude map, or don't
       if (trim(dpar%fg_filename(component)) == 'none') then
          constructor%amplitude = 0.d0
       else
          call read_bintab(trim(dpar%datadir)//trim(dpar%fg_filename(component)),&
               constructor%amplitude, npix, nmaps, nullval, anynull, header=header)
       end if

       ! Initialize pol_flag arrays
       constructor%nflag = 0
       constructor%pol_flag = 0

       do i = 1, constructor%nindices

          ! Load the poltype into the flag buffer, and store bit flags for each component
          flag_buffer = return_poltype_flag(dpar%fg_spec_poltype(component,i))
          constructor%nflag = size(flag_buffer)
          do j = 1, size(flag_buffer)
             constructor%pol_flag(i,j) = flag_buffer(j)
          end do

          ! Do we sample this index?
          constructor%sample_index(i)  = dpar%fg_samp_spec(component,i)

          ! Sample full sky or per-pixel?
          if (trim(dpar%fg_ind_region(component,i)) == 'fullsky') then
             constructor%index_mode = 1
          else if (trim(dpar%fg_ind_region(component,i)) == 'per-pixel') then
             constructor%index_mode = 2
          end if

          ! Define the lnl evaluation for each index
          constructor%lnl_type(i)      = dpar%fg_ind_lnl(component,i)

          ! Define MH step siz
          constructor%step_size(i) = dpar%fg_spec_step(component,i)

          ! Define prior for likelihood evaluation
          constructor%prior_type(i)    = dpar%fg_prior_type(component,i)
          constructor%gauss_prior(i,1) = dpar%fg_gauss(component,i,1)
          constructor%gauss_prior(i,2) = dpar%fg_gauss(component,i,2)
          constructor%uni_prior(i,1)   = dpar%fg_uni(component,i,1)
          constructor%uni_prior(i,2)   = dpar%fg_uni(component,i,2)

          ! Initialize spectral index maps, or don't
          if (trim(dpar%fg_spec_file(component,i)) == 'none') then
             constructor%indices(:,:,i)   = dpar%fg_init(component,i)
          else
             call read_bintab(trim(dpar%datadir)//trim(dpar%fg_spec_file(component,i)),&
                  constructor%indices(:,:,i), npix, nmaps, nullval, anynull, header=header)
          end if
       end do
    else if (trim(constructor%type) == 'cmb') then
       allocate(constructor%amplitude(0:npix-1,nmaps))
       constructor%nindices = 0

       ! Allocate general pol_type flag array
       allocate(constructor%nflag(1))
       allocate(constructor%pol_flag(1,3)) ! The three is here because that's the max number of poltypes we can handle

       ! Initialize pol_flag arrays
       constructor%nflag = 0
       constructor%pol_flag = 0

       ! Load the poltype into the flag buffer, and store bit flags for each component
       flag_buffer = return_poltype_flag(dpar%fg_spec_poltype(component,1))
       constructor%nflag = size(flag_buffer)
       do j = 1, size(flag_buffer)
          constructor%pol_flag(1,j) = flag_buffer(j)
       end do


       ! Initialize an amplitude map, or don't
       if (trim(dpar%fg_filename(component)) == 'none') then
          constructor%amplitude = 0.d0
       else
          call read_bintab(trim(dpar%datadir)//trim(dpar%fg_filename(component)),&
               constructor%amplitude, npix, nmaps, nullval, anynull, header=header)
       end if

    else if (trim(constructor%type) == 'template') then
       constructor%nindices = 0
      
       allocate(constructor%corr(nbands))
       allocate(constructor%amplitude(0:npix-1,nmaps)) ! Amplitude is given by the product of the template at the 
                                                       ! fit template amplitude (template_amplitudes)
       allocate(constructor%template_amplitudes(nbands,nmaps))
       allocate(constructor%template(0:npix-1,nmaps))

       ! Allocate general pol_type flag array
       allocate(constructor%nflag(1))
       allocate(constructor%pol_flag(1,3)) ! The three is here because that's the max number of poltypes we can handle

       ! Initialize pol_flag arrays
       constructor%nflag = 0
       constructor%pol_flag = 0

       ! Load the poltype into the flag buffer, and store bit flags for each component
       flag_buffer = return_poltype_flag(dpar%fg_spec_poltype(component,1))
       constructor%nflag = size(flag_buffer)
       do j = 1, size(flag_buffer)
          constructor%pol_flag(1,j) = flag_buffer(j)
       end do

       constructor%nfit                = dpar%fg_nfit(component) ! Also currently hardcoded
       constructor%corr                = dpar%fg_temp_corr(component,:)
       constructor%template_amplitudes = 0.d0

       ! Some error handling
       if (trim(dpar%fg_filename(component)) == 'none') then
          write(*,*) "Error: template filename == 'none' "
          stop
       else
          call read_bintab(trim(dpar%datadir)//trim(dpar%fg_filename(component)),&
               constructor%template, npix, nmaps, nullval, anynull, header=header)
       end if

    else if (trim(constructor%type) == 'hi_fit') then
       !------------------------------------------------------|
       ! Unique component!                                    |
       !------------------------------------------------------|
       ! signal model is given as follows                     |
       !                                                      |
       ! s_{nu,p} = A_{nu}*HI_{p}*B_{nu}(T_{d,p})             |
       !                                                      |
       ! So we have two things that we need to sample         |
       ! a template amplitude, and a spectral index per pixel |
       !------------------------------------------------------|
       constructor%nindices = 1
       
       allocate(constructor%corr(nbands))
       allocate(constructor%gauss_prior(1,2))
       allocate(constructor%uni_prior(1,2))
       allocate(constructor%sample_index(1))
       allocate(constructor%prior_type(1))
       allocate(constructor%lnl_type(1))
       allocate(constructor%index_mode(1))
       allocate(constructor%step_size(1))
       
       ! Allocate maps for the components
       allocate(constructor%template_amplitudes(nbands,nmaps))
       allocate(constructor%template(0:npix-1,nmaps))
       allocate(constructor%amplitude(0:npix-1,nmaps))
       allocate(constructor%indices(0:npix-1,nmaps,1))
       allocate(constructor%ind_label(constructor%nindices))

       constructor%ind_label = ['T']
       
       ! Reference frequency
       constructor%nu_ref              = dpar%fg_nu_ref(component) 

       constructor%nfit                = dpar%fg_nfit(component) ! Also currently hardcoded
       constructor%corr                = dpar%fg_temp_corr(component,:)
       constructor%template_amplitudes = 0.d0

       !HARD CODED
       constructor%step_size(1) = dpar%fg_spec_step(component,1)
       !!!!

       ! Allocate general pol_type flag array
       allocate(constructor%nflag(1))
       allocate(constructor%pol_flag(1,3)) ! The three is here because that's the max number of poltypes we can handle

       ! Initialize pol_flag arrays
       constructor%nflag = 0
       constructor%pol_flag = 0

       do i = 1, constructor%nindices

          ! Load the poltype into the flag buffer, and store bit flags for each component
          flag_buffer = return_poltype_flag(dpar%fg_spec_poltype(component,i))
          constructor%nflag = size(flag_buffer)
          do j = 1, size(flag_buffer)
             constructor%pol_flag(i,j) = flag_buffer(j)
          end do

          ! Define the lnl evaluation for each index
          constructor%lnl_type(i)      = dpar%fg_ind_lnl(component,i)

          ! Do we sample this index?
          constructor%sample_index(i)  = dpar%fg_samp_spec(component,i)

          ! Sample full sky or per-pixel?
          if (trim(dpar%fg_ind_region(component,i)) == 'fullsky') then
             constructor%index_mode = 1
          else if (trim(dpar%fg_ind_region(component,i)) == 'per-pixel') then
             constructor%index_mode = 2
          end if

          ! Define prior for likelihood evaluation
          constructor%prior_type(i)    = dpar%fg_prior_type(component,i)
          constructor%gauss_prior(i,1) = dpar%fg_gauss(component,i,1)
          constructor%gauss_prior(i,2) = dpar%fg_gauss(component,i,2)
          constructor%uni_prior(i,1)   = dpar%fg_uni(component,i,1)
          constructor%uni_prior(i,2)   = dpar%fg_uni(component,i,2)

          ! Initialize spectral index maps, or don't
          if (trim(dpar%fg_spec_file(component,i)) == 'none') then
             constructor%indices(:,:,i)   = dpar%fg_init(component,i)
          else
             write(*,*) component, i
             write(*,*) trim(dpar%fg_spec_file(component,i))
             call read_bintab(trim(dpar%datadir)//trim(dpar%fg_spec_file(component,i)),&
                  constructor%indices(:,:,i), npix, nmaps, nullval, anynull, header=header)
          end if
       end do

       ! Some error handling
       if (trim(dpar%fg_filename(component)) == 'none') then
          write(*,*) "Error: template filename == 'none' "
          stop
       else
          call read_bintab(trim(dpar%datadir)//trim(dpar%fg_filename(component)),&
               constructor%template, npix, nmaps, nullval, anynull, header=header)
       end if

    ! Some error handling
    else
       write(*,*) "Warning - unrecognized component type detected"
       stop
    end if
  end function constructor

  subroutine initialize_components(dpar)
    implicit none
    type(dang_params)             :: dpar
    integer(i4b)                  :: i

    write(*,*) 'Initializing components'

    allocate(component_list(dpar%ncomp))

    do i = 1, dpar%ncomp
       ! write(*,*) 'Initialize component ', i
       component_list(i)%p => dang_comps(dpar,i)
    end do

  end subroutine initialize_components

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

  function B_nu(nu,T)
    implicit none
    real(dp),       intent(in) :: nu, T
    real(dp)                   :: B_nu
    integer(i4b)               :: i
    ! Output in units of [W sr^-1 m^-2 Hz^-1]
    B_nu  = ((2.d0*h*(nu)**3.d0)/(c**2.d0))*(1.d0/(exp((h*nu)/(k_B*T))-1))
  end function B_nu

  function eval_signal(self, band, pix, map_n, index)
    implicit none
    class(dang_comps)                  :: self
    integer(i4b),           intent(in) :: band
    integer(i4b),           optional   :: pix
    integer(i4b),           optional   :: map_n
    integer(i4b)                       :: i
    real(dp), dimension(:), optional   :: index
    real(dp)                           :: eval_signal

    if (trim(self%type) == 'hi_fit') then
       eval_signal = self%template_amplitudes(band,map_n)*self%eval_sed(band,pix,map_n,index)
    else if (trim(self%type) == 'template') then
       eval_signal = self%template_amplitudes(band,map_n)*self%template(pix,map_n)
    else 
       eval_signal = self%amplitude(pix,map_n)*self%eval_sed(band,pix,map_n,index)
    end if

  end function eval_signal

  function eval_sed(self, band, pix, map_n, index)
    ! always computed in RJ units
    
    implicit none
    class(dang_comps)                  :: self
    integer(i4b),           intent(in) :: band
    integer(i4b),           optional   :: pix
    integer(i4b),           optional   :: map_n
    integer(i4b)                       :: i
    real(dp), dimension(:), optional   :: index
    real(dp)                           :: z, spectrum
    real(dp)                           :: eval_sed

    if (trim(self%type) == 'power-law') then
       spectrum = evaluate_powerlaw(self, band, pix, map_n, index)
    else if (trim(self%type) == 'mbb') then
       spectrum = evaluate_mbb(self, band, pix, map_n, index)
    else if (trim(self%type) == 'cmb') then
       spectrum = 1.0/a2t(bp(band))
    else if (trim(self%type) == 'template') then
       spectrum = 1.d0
    else if (trim(self%type) == 'hi_fit') then
       ! This isn't even in the right units!
       ! if (present(index)) then
       !    spectrum = self%template(pix,map_n)*planck(bp(band),index(1))
       ! else
       !    spectrum = self%template(pix,map_n)*planck(bp(band),self%indices(pix,map_n,1))
       ! end if
       spectrum = self%template(pix,map_n)*evaluate_hi_fit(self, band, pix, map_n, index)
    end if
    
    eval_sed = spectrum

  end function eval_sed

  function evaluate_hi_fit(self, band, pix, map_n, index)
    ! always computed in RJ units
    
    implicit none
    class(dang_comps)                  :: self
    integer(i4b),           intent(in) :: band
    integer(i4b),           optional   :: pix
    integer(i4b),           optional   :: map_n
    integer(i4b)                       :: i
    real(dp), dimension(:), optional   :: index
    real(dp)                           :: z, spectrum
    real(dp)                           :: evaluate_hi_fit

    spectrum = 0.d0

    if (bp(band)%id == 'delta') then
       if (present(index)) then
          spectrum = B_nu(bp(band)%nu_c,index(1))/compute_bnu_prime_RJ(bp(band)%nu_c)
       else
          spectrum = B_nu(bp(band)%nu_c,self%indices(pix,map_n,1))/compute_bnu_prime_RJ(bp(band)%nu_c)
       end if
    else
       if (present(index)) then
          do i = 1, bp(band)%n
             spectrum = spectrum + bp(band)%tau0(i)*B_nu(bp(band)%nu0(i),index(1))/&
                  & compute_bnu_prime_RJ(bp(band)%nu0(i))
          end do
       else
          do i = 1, bp(band)%n
             spectrum = spectrum + bp(band)%tau0(i)*B_nu(bp(band)%nu0(i),&
                  & self%indices(pix,map_n,1))/compute_bnu_prime_RJ(bp(band)%nu0(i))
          end do
       end if
    end if

    ! Make sure we move it to uK_RJ
    evaluate_hi_fit = spectrum*1e6

  end function evaluate_hi_fit

  function evaluate_powerlaw(self, band, pix, map_n, index)
    ! always computed in RJ units
    
    implicit none
    class(dang_comps)                  :: self
    integer(i4b),           intent(in) :: band
    integer(i4b),           optional   :: pix
    integer(i4b),           optional   :: map_n
    integer(i4b)                       :: i
    real(dp), dimension(:), optional   :: index
    real(dp)                           :: z, spectrum
    real(dp)                           :: evaluate_powerlaw

    spectrum = 0.d0

    if (bp(band)%id == 'delta') then
       if (present(index)) then
          spectrum = (bp(band)%nu_c/self%nu_ref)**index(1)
       else 
          spectrum = (bp(band)%nu_c/self%nu_ref)**self%indices(pix,map_n,1)
       end if
    else
       ! Compute for LFI bandpass
       if (present(index)) then
          do i = 1, bp(band)%n
             spectrum = spectrum + bp(band)%tau0(i)*(bp(band)%nu0(i)/self%nu_ref)**index(1)
          end do
       else 
          do i = 1, bp(band)%n
             spectrum = spectrum + bp(band)%tau0(i)*(bp(band)%nu0(i)/self%nu_ref)**self%indices(pix,map_n,1)
          end do
       end if
    end if

    evaluate_powerlaw = spectrum
    
  end function evaluate_powerlaw

  function evaluate_mbb(self, band, pix, map_n, index)
    ! always computed in RJ units
    
    implicit none
    class(dang_comps)                  :: self
    integer(i4b),           intent(in) :: band
    integer(i4b),           optional   :: pix
    integer(i4b),           optional   :: map_n
    integer(i4b)                       :: i
    real(dp), dimension(:), optional   :: index
    real(dp)                           :: z, spectrum
    real(dp)                           :: evaluate_mbb

    spectrum = 0.d0

    if (bp(band)%id == 'delta') then
       if (present(index)) then
          z = h / (k_B*index(2))
          spectrum = (exp(z*self%nu_ref)-1.d0) / &
               & (exp(z*bp(band)%nu_c)-1.d0) * (bp(band)%nu_c/(self%nu_ref))**(index(1)+1.d0)
       else
          z = h / (k_B*self%indices(pix,map_n,2))
          spectrum = (exp(z*self%nu_ref)-1.d0) / &
               & (exp(z*bp(band)%nu_c)-1.d0) * (bp(band)%nu_c/(self%nu_ref))**(self%indices(pix,map_n,1)+1.d0)
       end if
    else
       if (present(index)) then
          z = h / (k_B*index(2))
          do i = 1, bp(band)%n
             spectrum = spectrum + bp(band)%tau0(i)*(exp(z*self%nu_ref)-1.d0) / &
                  (exp(z*bp(band)%nu0(i))-1.d0) * (bp(band)%nu0(i)/(self%nu_ref))**(index(1)+1.d0)
          end do
       else
          z = h / (k_B*self%indices(pix,map_n,2))
          do i = 1, bp(band)%n
             spectrum = spectrum + bp(band)%tau0(i)*(exp(z*self%nu_ref)-1.d0) / &
                  (exp(z*bp(band)%nu0(i))-1.d0) * (bp(band)%nu0(i)/(self%nu_ref))**(self%indices(pix,map_n,1)+1.d0)
          end do
       end if
    end if

    evaluate_mbb = spectrum

  end function evaluate_mbb
    
end module dang_component_mod
