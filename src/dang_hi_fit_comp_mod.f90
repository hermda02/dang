module dang_hi_fit_comp_mod
  use dang_component_mod
  use dang_mixmat_1d_mod
  implicit none
  
  private
  public dang_hi_fit_comp

  type, extends (dang_comp) :: dang_hi_fit_comp
   contains
     procedure :: S => eval
  end type dang_hi_fit_comp

  interface dang_hi_fit_comp
     procedure constructor
  end interface dang_hi_fit_comp

contains

  function constructor(dpar,component)
    implicit none
    type(dang_params), intent(in) :: dpar
    integer(i4b),      intent(in) :: component
    class(dang_hi_fit_comp), pointer :: constructor

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
    mask_hi = .true.
    constructor%nindices = 1

    allocate(constructor%corr(nbands))
    allocate(constructor%gauss_prior(1,2))
    allocate(constructor%uni_prior(1,2))
    allocate(constructor%sample_index(1))
    allocate(constructor%sample_nside(1))
    allocate(constructor%prior_type(1))
    allocate(constructor%lnl_type(1))
    allocate(constructor%index_mode(1))

    allocate(constructor%tuned(constructor%nindices))
    allocate(constructor%step_size(1))

    ! Allocate maps for the components
    allocate(constructor%template_amplitudes(nbands,nmaps))
    allocate(constructor%template(0:npix-1,nmaps))
    allocate(constructor%temp_norm(nmaps))
    allocate(constructor%amplitude(0:npix-1,nmaps))
    allocate(constructor%indices(0:npix-1,nmaps,1))
    allocate(constructor%ind_label(constructor%nindices))

    constructor%ind_label = ['T']

    ! Reference frequency
    constructor%nu_ref              = dpar%fg_nu_ref(component) 

    constructor%nfit                = dpar%fg_nfit(component) ! Also currently hardcoded
    constructor%corr                = dpar%fg_temp_corr(component,:)
    constructor%template_amplitudes = 0.d0

    constructor%amplitude_file = trim(dpar%temp_amps(component))

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


       ! Define MH step size
       constructor%step_size(i) = dpar%fg_spec_step(component,1)
       constructor%tuned(i)     = .not. dpar%fg_spec_tune(component,i)

       ! Do we sample this index?
       constructor%sample_index(i)  = dpar%fg_samp_spec(component,i)

       ! What NSDIE?
       constructor%sample_nside(i)  = dpar%fg_samp_nside(component,i)

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
    class(dang_hi_fit_comp), intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pixel
    integer(i4b),            intent(in)           :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    integer(i4b)                                  :: i
    real(dp)                                      :: td, spectrum
    real(dp)                                      :: eval

    spectrum = 0.d0

    td = theta(1)

    if (present(nu)) then
       spectrum = B_nu(nu,td)/compute_bnu_prime_RJ(nu)
    else if (present(band)) then
       if (bp(band)%id == 'delta') then
          spectrum = B_nu(bp(band)%nu_c,td)/compute_bnu_prime_RJ(bp(band)%nu_c)
       else
          do i = 1, bp(band)%n
             if (bp(band)%nu0(i) == 0.d0) cycle 
             spectrum = spectrum + bp(band)%tau0(i)*B_nu(bp(band)%nu0(i),td)/&
                  & compute_bnu_prime_RJ(bp(band)%nu0(i))
          end do
       end if
    end if

    ! Make sure we move it to uK_RJ
    eval = spectrum*1e6

  end function eval
end module dang_hi_fit_comp_mod
