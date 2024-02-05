module dang_template_comp_mod
  use dang_component_mod
  use dang_mixmat_1d_mod
  implicit none
  
  private
  public dang_template_comp

  type, extends (dang_comp) :: dang_template_comp
   contains
     procedure :: S  => eval
  end type dang_template_comp

  interface dang_template_comp
     procedure constructor
  end interface dang_template_comp

contains

  function constructor(dpar,component)
    implicit none
    type(dang_params),      intent(in) :: dpar
    integer(i4b),           intent(in) :: component
    class(dang_template_comp), pointer :: constructor

    integer(i4b)                  :: i, j, k, count

    integer(i4b), allocatable, dimension(:) :: flag_buffer

    allocate(constructor)

    constructor%label            = trim(dpar%fg_label(component))
    constructor%type             = dpar%fg_type(component)
    constructor%cg_group         = dpar%fg_cg_group(component)

    if (trim(constructor%type) /= 'template' .and. trim(constructor%type) /= 'hi_fit' &
         & .and. trim(constructor%type) /= 'monopole') then
       constructor%sample_amplitude = dpar%fg_amp_samp(component)
    else if (trim(constructor%type) == 'T_cmb') then
       constructor%sample_amplitude = .false. ! Treat the CMB monopole as a uniform template
    else
       constructor%sample_amplitude = .true.
    end if

    write(*,*) 'Component label = '//trim(constructor%label)

    constructor%nindices = 0

    allocate(constructor%corr(nbands))
    allocate(constructor%amplitude(0:npix-1,nmaps)) ! Amplitude is given by the product of the template at the 
    ! fit template amplitude (template_amplitudes)
    allocate(constructor%template_amplitudes(nbands,nmaps))
    allocate(constructor%template(0:npix-1,nmaps))
    allocate(constructor%temp_norm(nmaps))

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

    do k = 1, nmaps
       constructor%temp_norm(k) = maxval(constructor%template(:,k))
       constructor%template(:,k) = constructor%template(:,k)/constructor%temp_norm(k)
    end do

  end function constructor
    
  function eval(self, nu, band, pol, theta, pixel)
    ! always computed in RJ units
    
    implicit none
    class(dang_template_comp), intent(in)           :: self
    real(dp),                  intent(in), optional :: nu
    integer(i4b),              intent(in), optional :: band
    integer(i4b),              intent(in), optional :: pixel
    integer(i4b),              intent(in)           :: pol
    real(dp), dimension(1:),   intent(in), optional :: theta
    integer(i4b)                       :: i
    real(dp)                           :: T, spectrum
    real(dp)                           :: eval

    if (self%corr(band)) then
       eval = self%template(pixel,pol)
    else
       eval = 0.d0
    end if

  end function eval
end module dang_template_comp_mod
