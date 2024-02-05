module dang_component_mod
  use healpix_types
  use pix_tools
  use fitstools
  use dang_util_mod
  use dang_param_mod
  use dang_bp_mod
  use dang_mixmat_mod
  implicit none
  
  public dang_comp, comp_ptr, component_list

  type, abstract :: dang_comp
     
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
     logical(lgt),      allocatable, dimension(:)     :: tuned            ! Do we tune this spectral index?
     integer(i4b),      allocatable, dimension(:)     :: sample_nside     ! At which nside do we sample?
     real(dp),          allocatable, dimension(:)     :: step_size        ! M-H step size

     ! Sky/template information
     real(dp),          allocatable, dimension(:,:)   :: amplitude           ! Amplitude maps
     real(dp),          allocatable, dimension(:,:)   :: template            ! Template map
     real(dp),          allocatable, dimension(:)     :: temp_norm           ! Template normalization factor
     real(dp),          allocatable, dimension(:,:)   :: template_amplitudes ! Template amplitudes

     integer(i4b),      allocatable, dimension(:)     :: index_mode          ! Fullsky/per-pixel
     real(dp),          allocatable, dimension(:,:,:) :: indices             ! Indices maps

     character(len=256)                               :: amplitude_file   ! File that holds the template amplitudes
     
     ! Prior information
     character(len=16), allocatable, dimension(:)     :: lnl_type
     character(len=16), allocatable, dimension(:)     :: prior_type
     real(dp),          allocatable, dimension(:,:)   :: gauss_prior
     real(dp),          allocatable, dimension(:,:)   :: uni_prior

     ! Poltype information
     integer(i4b),      allocatable, dimension(:)     :: nflag
     integer(i4b),      allocatable, dimension(:,:)   :: pol_flag

     type(mixmat_ptr), allocatable, dimension(:,:)    :: mixmat
     
   contains

     procedure                    :: evalSignal
     procedure(evalSED), deferred :: S
     
  end type dang_comp
  
  abstract interface
     function evalSED(self, nu, band, pol, theta, pixel)
       import dang_comp, dp, i4b
       class(dang_comp),        intent(in)           :: self
       real(dp),                intent(in), optional :: nu
       integer(i4b),            intent(in), optional :: band
       integer(i4b),            intent(in), optional :: pixel
       integer(i4b),            intent(in)           :: pol
       real(dp), dimension(1:), intent(in), optional :: theta
       real(dp)                                      :: evalSED
     end function evalSED
  end interface

  type comp_ptr
     class(dang_comp), pointer :: p => null()
  end type comp_ptr
  
  type(comp_ptr), allocatable, dimension(:) :: component_list 
  
contains 

  function evalSignal(self,band,pixel,pol,theta)
    implicit none
    class(dang_comp)                              :: self
    integer(i4b), intent(in)                      :: band, pixel, pol
    real(dp), dimension(1:), intent(in), optional :: theta
    real(dp)                                      :: evalSignal

    if (present(theta)) then
       evalSignal = self%amplitude(pixel,pol)*self%mixmat(band,pol)%p%eval(theta)
    else
       if ((trim(self%type) == 'monopole') .or. (trim(self%type) == 'template')) then
          evalSignal = self%template_amplitudes(band,pol)*self%template(pixel,pol)
       else
          evalSignal = self%amplitude(pixel,pol)*self%mixmat(band,pol)%p%eval(self%indices(pixel,pol,:))
       end if
    end if
       
  end function evalSignal
  
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
         if (bp%nu0(i) == 0.d0) cycle 
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

end module dang_component_mod
