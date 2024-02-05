module dang_cmb_comp_mod
  use dang_component_mod
  use dang_mixmat_1d_mod
  implicit none
  
  private
  public dang_cmb_comp

  type, extends (dang_comp) :: dang_cmb_comp
   contains
     procedure :: S  => eval
  end type dang_cmb_comp

  interface dang_cmb_comp
     procedure constructor
  end interface dang_cmb_comp

contains

  function constructor(dpar,component)
    implicit none
    type(dang_params), intent(in) :: dpar
    integer(i4b),      intent(in) :: component
    class(dang_cmb_comp), pointer :: constructor
    
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
    class(dang_cmb_comp),    intent(in)           :: self
    real(dp),                intent(in), optional :: nu
    integer(i4b),            intent(in), optional :: band
    integer(i4b),            intent(in), optional :: pixel
    integer(i4b),            intent(in)           :: pol
    real(dp), dimension(1:), intent(in), optional :: theta
    integer(i4b)                       :: i
    real(dp)                           :: T, spectrum
    real(dp)                           :: eval

    spectrum = 0.d0
    
    if (present(nu)) then
       spectrum = 1.d0/a2t_nu(nu)
    else if (present(band)) then
       if (bp(band)%id == 'delta') then
          spectrum = 1.d0/a2t_nu(bp(band)%nu_c)
       else
          do i = 1, bp(band)%n
             if (bp(band)%nu0(i) == 0.d0) cycle 
             spectrum = spectrum + bp(band)%tau0(i)/a2t_nu(bp(band)%nu0(i))
          end do
       end if
    end if

    eval = spectrum*1e6
    
  end function eval
end module dang_cmb_comp_mod
