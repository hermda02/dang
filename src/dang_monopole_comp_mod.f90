module dang_monopole_comp_mod
  use dang_component_mod
  use dang_mixmat_1d_mod
  implicit none
  
  private
  public dang_monopole_comp

  type, extends (dang_comp) :: dang_monopole_comp
   contains
     procedure :: S  => eval
  end type dang_monopole_comp

  interface dang_monopole_comp
     procedure constructor
  end interface dang_monopole_comp

contains

  function constructor(dpar,component)
    implicit none
    type(dang_params),      intent(in) :: dpar
    integer(i4b),           intent(in) :: component
    class(dang_monopole_comp), pointer :: constructor

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

    ! Allocate general pol_type flag array
    allocate(constructor%nflag(1))
    allocate(constructor%pol_flag(1,3)) ! The three is here because that's the max number of poltypes we can handle

    constructor%nindices = 0

    allocate(constructor%corr(nbands))
    allocate(constructor%amplitude(0:npix-1,nmaps)) ! Amplitude is given by the product of the template at the 
    ! fit template amplitude (template_amplitudes)
    allocate(constructor%template_amplitudes(nbands,nmaps))
    allocate(constructor%template(0:npix-1,nmaps))

    constructor%nfit                = dpar%fg_nfit(component) ! Also currently hardcoded
    constructor%corr                = dpar%fg_temp_corr(component,:)
    constructor%template_amplitudes = 0.d0
    ! Set polarized intensity offset map to 0
    constructor%template            = 0.d0
    ! Set total intensity offset map to 1
    constructor%template(:,1)       = 1.d0

  end function constructor
    
  function eval(self, nu, band, pol, theta, pixel)
    ! always computed in RJ units
    
    implicit none
    class(dang_monopole_comp), intent(in)           :: self
    real(dp),                  intent(in), optional :: nu
    integer(i4b),              intent(in), optional :: band
    integer(i4b),              intent(in), optional :: pixel
    integer(i4b),              intent(in)           :: pol
    real(dp), dimension(1:),   intent(in), optional :: theta
    integer(i4b)                       :: i
    real(dp)                           :: T, spectrum
    real(dp)                           :: eval

    if (self%corr(band)) then
       eval = 1.d0
    else
       eval = 0.d0
    end if

  end function eval
end module dang_monopole_comp_mod
