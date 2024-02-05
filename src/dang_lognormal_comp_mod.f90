module dang_lognormal_comp_mod
  use dang_component_mod
  use dang_mixmat_1d_mod
  implicit none
  
  private
  public dang_lognormal_comp


  type, extends (dang_comp) :: dang_lognormal_comp
   contains
     procedure :: S  => eval
  end type dang_lognormal_comp

  interface dang_lognormal_comp
     procedure constructor
  end interface dang_lognormal_comp

contains

  function constructor(dpar,component)
    implicit none
    type(dang_params), intent(in) :: dpar
    integer(i4b),      intent(in) :: component
    class(dang_lognormal_comp), pointer    :: constructor

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
    

       ! Allocate arrays to appropriate size for each component type
       constructor%nindices = 2
       allocate(constructor%gauss_prior(constructor%nindices,2))
       allocate(constructor%uni_prior(constructor%nindices,2))
       allocate(constructor%sample_index(constructor%nindices))
       allocate(constructor%sample_nside(constructor%nindices))
       allocate(constructor%prior_type(constructor%nindices))
       allocate(constructor%lnl_type(constructor%nindices))
       allocate(constructor%index_mode(constructor%nindices))

       allocate(constructor%step_size(constructor%nindices))
       allocate(constructor%tuned(constructor%nindices))

       ! Allocate general pol_type flag array
       allocate(constructor%nflag(constructor%nindices))
       allocate(constructor%pol_flag(constructor%nindices,3)) ! The three is here because that's the max number of poltypes we can handle

       ! Allocate maps for the components
       allocate(constructor%amplitude(0:npix-1,nmaps))
       allocate(constructor%indices(0:npix-1,nmaps,constructor%nindices))

       allocate(constructor%ind_label(constructor%nindices))

       constructor%ind_label = ['nu_p ', 'w_ame']

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

          ! What NSDIE?
          constructor%sample_nside(i)  = dpar%fg_samp_nside(component,i)

          ! Sample full sky or per-pixel?
          if (trim(dpar%fg_ind_region(component,i)) == 'fullsky') then
             constructor%index_mode(i) = 1
          else if (trim(dpar%fg_ind_region(component,i)) == 'per-pixel') then
             constructor%index_mode(i) = 2
          end if

          ! Define the lnl evaluation for each index
          constructor%lnl_type(i)      = dpar%fg_ind_lnl(component,i)

          ! Define MH step size
          constructor%step_size(i) = dpar%fg_spec_step(component,i)
          constructor%tuned(i)     = .not. dpar%fg_spec_tune(component,i)

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
    class(dang_lognormal_comp), intent(in)        :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pixel
    integer(i4b),            intent(in)           :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    integer(i4b)                       :: i
    real(dp)                           :: nu_p, w_ame, spectrum
    real(dp)                           :: eval

    spectrum = 0.d0

    ! Extract indices
    nu_p  = theta(1)
    w_ame = theta(2)
    
    if (present(nu)) then
       spectrum = exp(-0.5*(log(nu/(nu_p*1e9))/w_ame)**2)*(self%nu_ref/nu)**2
    else if (present(band)) then
       ! Evaluate
       if (bp(band)%id == 'delta') then
          spectrum = exp(-0.5*(log(bp(band)%nu_c/(nu_p*1e9))/w_ame)**2)*(self%nu_ref/bp(band)%nu_c)**2
       else
          do i = 1, bp(band)%n
             if (bp(band)%nu0(i) == 0.d0) cycle 
             spectrum = spectrum + bp(band)%tau0(i)*&
                  & exp(-0.5*(log(bp(band)%nu0(i)/(nu_p*1e9))/w_ame)**2)*(self%nu_ref/bp(band)%nu0(i))**2
          end do
       end if
    end if

    eval = spectrum
    
  end function eval
end module dang_lognormal_comp_mod
