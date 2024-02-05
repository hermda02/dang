module dang_T_cmb_comp_mod
  use dang_component_mod
  use dang_mixmat_1d_mod
  implicit none

  private
  public dang_T_cmb_comp


  type, extends (dang_comp) :: dang_T_cmb_comp
   contains
     procedure :: S  => eval
  end type dang_T_cmb_comp

  interface dang_T_cmb_comp
     procedure constructor
  end interface dang_T_cmb_comp

contains

  function constructor(dpar,component)
    implicit none
    type(dang_params), intent(in) :: dpar
    integer(i4b),      intent(in) :: component
    class(dang_T_cmb_comp), pointer :: constructor
    
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

    constructor%nindices = 1

    allocate(constructor%corr(nbands))
    allocate(constructor%gauss_prior(1,2))
    allocate(constructor%uni_prior(1,2))
    allocate(constructor%sample_index(1))
    allocate(constructor%sample_nside(1))
    allocate(constructor%prior_type(1))
    allocate(constructor%lnl_type(1))
    allocate(constructor%index_mode(1))

    allocate(constructor%step_size(1))
    allocate(constructor%tuned(constructor%nindices))

    ! Allocate maps for the components
    allocate(constructor%amplitude(0:npix-1,nmaps))
    allocate(constructor%indices(0:npix-1,nmaps,1))
    allocate(constructor%ind_label(constructor%nindices))

    constructor%ind_label = ['T']

    ! Allocate general pol_type flag array
    allocate(constructor%nflag(1))
    allocate(constructor%pol_flag(1,3)) ! The three is here because that's the max number of poltypes we can handle

    ! Load the poltype into the flag buffer, and store bit flags for each component
    constructor%nflag = 1
    constructor%pol_flag(1,1) = 1 ! Only total intensity
    constructor%index_mode(1) = 1 ! Only assume the temperature is fullsky

    ! Treat the CMB like a template, essentially
    constructor%amplitude = 1.d0

    ! Do we sample this index?
    constructor%sample_index(1)  = dpar%fg_samp_spec(component,1)

    ! What NSDIE?
    constructor%sample_nside(1)  = dpar%fg_samp_nside(component,1)

    ! Define the lnl evaluation for each index
    constructor%lnl_type(1)      = dpar%fg_ind_lnl(component,1)

    ! Define MH step size
    constructor%step_size(1)     = dpar%fg_spec_step(component,1)
    constructor%tuned(1)         = .not. dpar%fg_spec_tune(component,i)

    ! Define prior for likelihood evaluation
    constructor%prior_type(1)    = dpar%fg_prior_type(component,1) 
    constructor%gauss_prior(1,1) = dpar%fg_gauss(component,1,1)
    constructor%gauss_prior(1,2) = dpar%fg_gauss(component,1,2)
    constructor%uni_prior(1,1)   = dpar%fg_uni(component,1,1)
    constructor%uni_prior(1,2)   = dpar%fg_uni(component,1,2)

    constructor%indices(:,:,1)   = dpar%fg_init(component,1)

    allocate(constructor%mixmat(nbands,3))
    do k = 1, 3
       do j = 1, nbands
          constructor%mixmat(j,k)%p => dang_mixmat_1d(constructor, bp(j), k)
       end do
    end do

  end function constructor

  function eval(self, nu, band, pol, theta, pixel)
    ! always computed in RJ units

    implicit none
    class(dang_T_cmb_comp),  intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pixel
    integer(i4b),            intent(in)           :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    integer(i4b)                       :: i
    real(dp)                           :: T, spectrum
    real(dp)                           :: eval

    spectrum = 0.d0

    T = theta(1)

    if (present(nu)) then
       spectrum = B_nu(nu,T)/compute_bnu_prime_RJ(nu)
    else if (present(band)) then
       if (bp(band)%id == 'delta') then
          spectrum = B_nu(bp(band)%nu_c,T)/compute_bnu_prime_RJ(bp(band)%nu_c)
       else
          do i = 1, bp(band)%n
             if (bp(band)%nu0(i) == 0.d0) cycle 
             spectrum = spectrum + bp(band)%tau0(i)*B_nu(bp(band)%nu0(i),T)/&
                  & compute_bnu_prime_RJ(bp(band)%nu0(i))
          end do
       end if
    end if

    eval = spectrum*1e6

  end function eval
end module dang_T_cmb_comp_mod
